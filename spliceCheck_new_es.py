import os
import requests
import json
# import urllib.parse
import sys
import operator
import time
import asyncio

nucleotide_change = {"A": "T", "T": "A", "C": "G", "G": "C"}

def create_mes_dict(wtseq, mutseq, prime, strand):
	""" Create a dictionary of sequences 10 bases downstream and upstream of the sequence """
	mes_dict = {}
	# print(wtseq, mutseq, len(wtseq), "\n")
	if prime == 5:
		for i in range(0, len(wtseq)):
			if len(wtseq[i:i+9]) != 9:
				break
			nine_mer_wt = wtseq[i:i+9]
			nine_mer_mut = mutseq[i:i+9]
			mes5_wt = mes5(nine_mer_wt)
			mes5_mut = mes5(nine_mer_mut)
			if mes5_wt < 0:
				mes5_wt = 0
			if mes5_mut < 0:
				mes5_mut = 0
			mes_dict[i - 11] = (mes5_wt, mes5_mut, round(mes5_mut - mes5_wt, 2))
			# print(i - 11)
			# print(nine_mer_wt, nine_mer_mut)
			# print(mes5_wt, mes5_mut)
			# print("--------------------------------------")
	elif prime == 3:
		for i in range(0, len(wtseq)):
			if len(wtseq[i:i+23]) != 23 or i - 11 > 11:
				break
			twentythree_mer_wt = wtseq[i:i+23]
			twentythree_mer_mut = mutseq[i:i+23]
			mes3_wt = mes3(twentythree_mer_wt)
			mes3_mut = mes3(twentythree_mer_mut)
			if mes3_wt < 0:
				mes3_wt = 0
			if mes3_mut < 0:
				mes3_mut = 0
			mes_dict[i - 11] = (mes3_wt, mes3_mut, round(mes3_mut - mes3_wt, 2))
			# print(i - 11)
			# print(twentythree_mer_wt, twentythree_mer_mut)
			# print(mes3_wt, mes3_mut)
			# print("--------------------------------------")
	return mes_dict

def mes5(seq):
	""" Call 5 prime end MaxEntScan algorithm http://genes.mit.edu/burgelab/maxent/download/"""
	# Seq example: CAGgtaagt
	if len(seq) != 9:
		raise Exception

	# Run maxentscan
	output = (os.popen('cd tools/maxentscan; perl score5_mod.pl %s' % (seq)).readline())

	# If the score is 0 it is outputted as an empty string
	if output == '':
		return 0.00
	else:
		return float(output)

def mes3(seq):
	""" Call 5 prime end MaxEntScan algorithm http://genes.mit.edu/burgelab/maxent/download/"""
	# Seq example: ttccaaacgaacttttgtagGGA
	if len(seq) != 23:
		raise Exception

	# Run maxentscan
	output = (os.popen('cd tools/maxentscan; perl score3_mod.pl %s' % (seq)).readline())

	# If the score is 0 it is outputted as an empty string
	if output == '':
		return 0.00
	else:
		return float(output)

def spliceai(gnomad):
	try:
		# Retrieve PubMed IDs
		server = "https://spliceailookup-api.broadinstitute.org/"
		get = "spliceai/"
		params = {"hg": 38, "distance":200, "variant":gnomad}	# change to 5000 for ATM project
		r = requests.get(server + get, params=params, headers={"Content-Type": "application/json"}, verify=False)

		if not r.ok:
			r.raise_for_status()
			return None
		res = r.json()
		return str(res["scores"][0])
	except Exception as e:
		print(e)
		return []

def determine_impacts(gen_chr, gen_start, hgvs, consequence, sift_score, polyphen_score, maxentscan_ref, maxentscan_alt, maxentscan_diff, spliceai_pred, mes_dict_5, mes_dict_3, strand):

	# data structures to help determine aso amenability
	mes_analyses = {}
	spliceai_analyses = {}

	# thresholds
	sift_threshold = 0.05
	polyphen_threshold = 0.85
	maxentscan_threshold = 4
	spliceai_threshold = 0.05

	# check the coding impact
	coding_impact = ""
	if consequence == "stop_gained" or "stop" in consequence:	# Stop gain mutation
		coding_impact = "VEP classified this variant as a stop-gain mutation."
	elif consequence == "splice_donor_variant" or consequence == "splice_acceptor_variant" or consequence == "splice_region_variant":	# splice region variant
		coding_impact = "VEP classified this variant as a {}.".format(consequence)
	else:	# missense
		score_pred = ""
		if type(sift_score) is list and type(polyphen_score) is list:	# no sift and polyphen
			score_pred = "No SIFT or PolyPhen score found."
		elif type(sift_score) is not list and type(polyphen_score) is not list:	# have both sift and polyphen
			if sift_score < sift_threshold and polyphen_score > polyphen_threshold:
				score_pred = "The scores indicate the mutation is predicted to be damaging/deleterious."
			elif sift_score > sift_threshold and polyphen_score < polyphen_threshold:
				score_pred = "The scores indicate that the mutation is predicted to be benign/tolerated."
			else:
				score_pred = "The scores indicate that the coding impact of this mutation is confusing."
		elif type(sift_score) is not list:
			if sift_score < sift_threshold:
				score_pred = "The SIFT score indicates that the mutation is predicted to be damaging/deleterious but without a PolyPhen score it deserves a further look."
			else:
				score_pred = "The SIFT score indicates that the mutation is predicted to be benign/tolerated but without a PolyPhen score it deserves a further look."
		else:
			if polyphen_score > polyphen_threshold:
				score_pred = "The PolyPhen score indicates that the mutation is predicted to be damaging/deleterious but without a SIFT score it deserves a further look."
			else:
				score_pred = "The PolyPhen score indicates that the mutation is predicted to be benign/tolerated but without a SIFT score it deserves a further look."
		coding_impact = "VEP classified this variant as a(n) {}. {}".format(consequence, score_pred)

	# check the splicing impact - MES/SpliceAI
	splicing_impact = ""
	maxentscan_res = ""
	spliceai = ""

	if consequence == "splice_donor_variant" or consequence == "splice_region_variant":
		dicts = [(mes_dict_5, "5")]
	elif consequence == "splice_acceptor_variant":
		dicts = [(mes_dict_3, "3")]
	else:
		dicts = [(mes_dict_5, "5"), (mes_dict_3, "3")]

	results = {}
	for dict_tuple in dicts:
		dict_to_use = dict_tuple[0]
		if dict_tuple[1] == "5":
			splice_qual, splice_qual_num = "donor", 1
		else:
			splice_qual, splice_qual_num = "acceptor", 0
		for key in dict_to_use:
			if dict_to_use[key][2] != 0:	# if the delta is not 0
				ref, alt, delta = dict_to_use[key][0], dict_to_use[key][1], dict_to_use[key][2]
				dict_key = "{}:{}:{}:{}".format(str(key), str(ref), str(alt), str(delta))

				# check native loss
				if delta < 0:
					if alt < 6.2:
						if delta <= (-1.15):
							results[dict_key] = "MES indicates a high probability (ref = {}, alt = {}, delta = {}) that a splice {} site is being weakened or destroyed.".format(ref, alt, delta, splice_qual)
							mes_analyses[dict_key] = (2, splice_qual_num, 0)
						else:
							results[dict_key] = "MES indicates a moderate probability (ref = {}, alt = {}, delta = {}) that a splice {} site is being weakened or destroyed.".format(ref, alt, delta, splice_qual)
							mes_analyses[dict_key] = (1, splice_qual_num, 0)
					elif alt >= 6.2 and alt <= 8.5:
						if delta <= (-1.15):
							results[dict_key] = "MES indicates a moderate probability (ref = {}, alt = {}, delta = {}) that a splice {} site is being weakened or destroyed.".format(ref, alt, delta, splice_qual)
							mes_analyses[dict_key] = (1, splice_qual_num, 0)
						else:
							results[dict_key] = "MES indicates a low probability (ref = {}, alt = {}, delta = {}) that a splice {} site is being weakened or destroyed.".format(ref, alt, delta, splice_qual)
							mes_analyses[dict_key] = (0, splice_qual_num, 0)
					else:
						results[dict_key] = "MES indicates a low probability (ref = {}, alt = {}, delta = {}) that a splice {} site is being weakened or destroyed.".format(ref, alt, delta, splice_qual)
						mes_analyses[dict_key] = (0, splice_qual_num, 0)
				# check de novo gain
				else:
					if alt > 8.5:
						results[dict_key] = "MES indicates a high probability (ref = {}, alt = {}, delta = {}) that a de novo/cryptic splice {} site is being strengthened or created.".format(ref, alt, delta, splice_qual)
						mes_analyses[dict_key] = (2, splice_qual_num, 1)
					elif alt < 6.2:
						results[dict_key] = "MES indicates a low probability (ref = {}, alt = {}, delta = {}) that a de novo/cryptic splice {} site is being strengthened or created.".format(ref, alt, delta, splice_qual)
						mes_analyses[dict_key] = (0, splice_qual_num, 1)
					else:
						if "+" in hgvs or "-" in hgvs:	# intronic
							for other_key in dict_to_use:	# compare current ref/alt/delta to the other splice sites
								if other_key < key:	# make sure we are only looking at upstream splice sites
									if alt > dict_to_use[other_key][1]:
										results[dict_key] = "MES indicates a moderate probability (ref = {}, alt = {}, delta = {}) that a de novo/cryptic splice {} site is being strengthened or created.".format(ref, alt, delta, splice_qual)
										mes_analyses[dict_key] = (1, splice_qual_num, 1)
									else:
										results[dict_key] = "MES indicates a low probability (ref = {}, alt = {}, delta = {}) that a de novo/cryptic splice {} site is being strengthened or created.".format(ref, alt, delta, splice_qual)
										mes_analyses[dict_key] = (0, splice_qual_num, 1)
						else:	# exonic
							for other_key in dict_to_use:	# compare current ref/alt/delta to the other splice sites
								if other_key > key:	# make sure we are only looking at downstream splice sites
									if alt > dict_to_use[other_key][1]:
										results[dict_key] = "MES indicates a moderate probability (ref = {}, alt = {}, delta = {}) that a de novo/cryptic splice {} site is being strengthened or created.".format(ref, alt, delta, splice_qual)
										mes_analyses[dict_key] = (1, splice_qual_num, 1)
									else:
										results[dict_key] = "MES indicates a low probability (ref = {}, alt = {}, delta = {}) that a de novo/cryptic splice {} site is being strengthened or created.".format(ref, alt, delta, splice_qual)
										mes_analyses[dict_key] = (0, splice_qual_num, 1)

	if len(results.keys()) == 0:	# nothing found from MES Perl scripts - try VEP
		if (type(maxentscan_ref) is not list) and (type(maxentscan_alt) is not list) and (type(maxentscan_diff) is not list):
			maxentscan_diff = maxentscan_diff * -1	# change the delta to match our calculation (alt - ref) instead of VEP's (ref - alt)
			dict_key = "{}:{}:{}".format(str(maxentscan_ref), str(maxentscan_alt), str(maxentscan_diff))
			if maxentscan_diff < 0:
				if maxentscan_ref < 6.2:
					if maxentscan_diff <= (-1.15):
						results[dict_key] = "MES indicates a high probability (ref = {}, alt = {}, delta = {}) that a splice site is being weakened or destroyed.".format(maxentscan_ref, maxentscan_alt, maxentscan_diff)
						mes_analyses[dict_key] = (2, 2, 0)
					else:
						results[dict_key] = "MES indicates a moderate probability (ref = {}, alt = {}, delta = {}) that a splice site is being weakened or destroyed.".format(maxentscan_ref, maxentscan_alt, maxentscan_diff)
						mes_analyses[dict_key] = (1, 2, 0)
				elif maxentscan_alt >= 6.2 and maxentscan_alt <= 8.5:
					if maxentscan_diff <= (-1.15):
						results[dict_key] = "MES indicates a moderate probability (ref = {}, alt = {}, delta = {}) that a splice site is being weakened or destroyed.".format(maxentscan_ref, maxentscan_alt, maxentscan_diff)
						mes_analyses[dict_key] = (1, 2, 0)
					else:
						results[dict_key] = "MES indicates a low probability (ref = {}, alt = {}, delta = {}) that a splice site is being weakened or destroyed.".format(maxentscan_ref, maxentscan_alt, maxentscan_diff)
						mes_analyses[dict_key] = (0, 2, 0)
				elif maxentscan_diff <= -4:
					results[dict_key] = "MES indicates a moderate probability (ref = {}, alt = {}, delta = {}) that a splice site is being weakened or destroyed.".format(maxentscan_ref, maxentscan_alt, maxentscan_diff)
					mes_analyses[dict_key] = (1, 2, 0)
				else:
					results[dict_key] = "MES indicates a low probability (ref = {}, alt = {}, delta = {}) that a splice site is being weakened or destroyed.".format(maxentscan_ref, maxentscan_alt, maxentscan_diff)
					mes_analyses[dict_key] = (0, 2, 0)
			else:
				if maxentscan_alt > 8.5:
					results[dict_key] = "MES indicates a high probability (ref = {}, alt = {}, delta = {}) that a de novo/cryptic splice site is being strengthened or created.".format(maxentscan_ref, maxentscan_alt, maxentscan_diff)
					mes_analyses[dict_key] = (2, 2, 1)
				elif maxentscan_alt < 6.2:
					results[dict_key] = "MES indicates a low probability (ref = {}, alt = {}, delta = {}) that a de novo/cryptic splice site is being strengthened or created.".format(maxentscan_ref, maxentscan_alt, maxentscan_diff)
					mes_analyses[dict_key] = (0, 2, 1)
				else:
					results[dict_key] = "MES indicates a low/moderate probability (ref = {}, alt = {}, delta = {}) that a de novo/cryptic splice site is being strengthened or created.".format(maxentscan_ref, maxentscan_alt, maxentscan_diff)
					mes_analyses[dict_key] = (1, 2, 1)

	maxentscan_res = ""
	gen_start = str(int(gen_start) + 1)
	for key in results:
		if len(key.split(":")) == 4:
			pos = key.split(":")[0]
			delta_pos = str(int(gen_start) + int(pos))
			# minor string manipulation
			sign = "+"
			if int(pos) < 0:
				sign = "-"
				if len(pos) >= 2:
					pos = pos[1:]
			maxentscan_res += "At position {}:{} (={} {} {}), {} ".format(gen_chr, delta_pos, gen_start, sign, pos, results[key])
		else:
			maxentscan_res += "From VEP results, {} ".format(results[key])
	if len(maxentscan_res) == 0:
		maxentscan_res = "No notable MaxEntScan results found."

	# spliceAI analysis
	if type(spliceai_pred) is list:
		spliceai = "No SpliceAI prediction found."
	else:
		if "---" in spliceai_pred:
			spliceai_info = spliceai_pred.strip().split("---")[-1].split("|")
		else:
			spliceai_info = spliceai_pred.strip().split("|")
		acceptor_gain, acceptor_loss, donor_gain, donor_loss = float(spliceai_info[1]), float(spliceai_info[2]), float(spliceai_info[3]), float(spliceai_info[4])
		acc_gain_pos, acc_loss_pos, donor_gain_pos, donor_loss_pos = strand * float(spliceai_info[5]), strand * float(spliceai_info[6]), strand * float(spliceai_info[7]), strand * float(spliceai_info[8])
		if acceptor_gain >= spliceai_threshold:
			spliceai = spliceai + "SpliceAI indicates the probability that the position {}bp away is used as a splice acceptor increases by {}. ".format(acc_gain_pos, acceptor_gain)
			spliceai_analyses["{}:{}".format(acc_gain_pos, acceptor_gain)] = (0, 1)
		if acceptor_loss >= spliceai_threshold:
			spliceai = spliceai + "SpliceAI indicates the probability that the position {}bp away is used as a splice acceptor decreases by {}. ".format(acc_loss_pos, acceptor_loss)
			spliceai_analyses["{}:{}".format(acc_loss_pos, acceptor_loss)] = (0, 0)
		if donor_gain >= spliceai_threshold:
			spliceai = spliceai + "SpliceAI indicates the probability that the position {}bp away is used as a splice donor increases by {}. ".format(donor_gain_pos, donor_gain)
			spliceai_analyses["{}:{}".format(donor_gain_pos, donor_gain)] = (1, 1)
		if donor_loss >= spliceai_threshold:
			spliceai = spliceai + "SpliceAI indicates the probability that the position {}bp away is used as a splice donor decreases by {}. ".format(donor_loss_pos, donor_loss)
			spliceai_analyses["{}:{}".format(donor_loss_pos, donor_loss)] = (1, 0)
		if len(spliceai) == 0:
			spliceai = "No notable results from SpliceAI."

	splicing_impact = "{} {}".format(maxentscan_res, spliceai)

	return coding_impact, splicing_impact, consequence, mes_analyses, spliceai_analyses

def mes5_runner(gen_start, gen_chr, mut, strand, hgvs, indel_length=0, indel_bp=""):
	# run 5' maxentscan
	mes_dict_5 = {}
	mes5lowdiff, mes5highdiff = 13, 19	# THESE VALUES WERE DETERMINED THROUGH TRIAL AND ERROR
	if strand == -1:
		mes5lowdiff, mes5highdiff = mes5highdiff - 1, mes5lowdiff + 1

	mes5highdiff += indel_length
	try:
		mes5low = str(int(gen_start) - mes5lowdiff)
		mes5high = str(int(gen_start) + mes5highdiff)
		# mes5range = "chr{}:{}-{}".format(gen_chr, mes5low, mes5high)
		try:
			server = "https://api.genome.ucsc.edu"
			get = "/getData/sequence"
			params = {"chrom": "chr{}".format(gen_chr), "genome": "hg38", "start": mes5low, "end": mes5high}
			r = requests.get(server + get, params=params, headers={"Content-Type": "application/json"}, verify=True)
			if not r.ok:
				r.raise_for_status()
				return {}
			decoded = r.json()
			res_sequence = decoded["dna"]
		except Exception as e:
			print("Could not retrieve sequence from UCSC: ", e)
			return {}
		if res_sequence:
			if strand == -1:
				new_res_sequence = ""
				for char in res_sequence[::-1].strip():
					new_res_sequence += nucleotide_change[char]
				seq9wt = new_res_sequence[0:31].lower()
				mes5lowdiff, mes5highdiff = mes5highdiff - 1, mes5lowdiff + 1
				res_sequence = new_res_sequence
			else:
				seq9wt = res_sequence[0:31].lower()

			if "del" in hgvs:
				seq9mut = "{}{}".format(res_sequence[0:mes5lowdiff], res_sequence[mes5lowdiff+indel_length:len(res_sequence)-1]).lower()
			elif "ins" in hgvs:
				seq9wt += res_sequence[31:31+indel_length].lower()
				seq9mut = "{}{}{}".format(res_sequence[0:mes5lowdiff], indel_bp, res_sequence[mes5lowdiff:mes5lowdiff+mes5highdiff-len(indel_bp)-1]).lower()
			elif "dup" in hgvs:
				seq9wt += res_sequence[31:31+indel_length].lower()
				seq9mut = "{}{}{}".format(res_sequence[0:mes5lowdiff], indel_bp, res_sequence[mes5lowdiff:mes5lowdiff+mes5highdiff-len(indel_bp)-1]).lower()
			else:
				seq9mut = "{}{}{}".format(seq9wt[0:mes5lowdiff], mut, seq9wt[mes5lowdiff+1:mes5lowdiff+mes5highdiff]).lower()
			# print(res_sequence.lower(), mes5lowdiff, mes5highdiff, indel_length, indel_bp)
			# print(seq9wt, len(seq9wt))
			# print(seq9mut, len(seq9mut), "\n")
			return create_mes_dict(seq9wt, seq9mut, 5, strand)
		else:
			print("Something went wrong with 5' maxentscan. Please re-enter your coordinates.")
			return {}
	except Exception as e:
		print("Could not run 5' maxentscan: ", e)
		wt5mes = [False]
		mut5mes = [False]
		return {}

def mes3_runner(gen_start, gen_chr, mut, strand, hgvs, indel_length=0, indel_bp=""):
	mes_dict_3 = {}
	mes23lowdiff, mes23highdiff = 31, 21
	if strand == -1:
		mes23lowdiff, mes23highdiff = mes23highdiff - 1, mes23lowdiff + 1

	mes23highdiff += indel_length
	try:
		mes23low = str(int(gen_start) - mes23lowdiff)
		mes23high = str(int(gen_start) + mes23highdiff)
		# mes23range = "chr{}:{}-{}".format(gen_chr, mes23low, mes23high)
		try:
			server = "https://api.genome.ucsc.edu"
			get = "/getData/sequence"
			params = {"chrom": "chr{}".format(gen_chr), "genome": "hg38", "start": mes23low, "end": mes23high}
			r = requests.get(server + get, params=params, headers={"Content-Type": "application/json"})
			if not r.ok:
				r.raise_for_status()
				return {}
			decoded = r.json()
			res_sequence = decoded["dna"]
		except Exception as e:
			print("Could not retrieve sequence from UCSC: ", e)
			return {}
		if res_sequence:
			if strand == -1:
				new_res_sequence = ""
				for char in res_sequence[::-1].strip():
					new_res_sequence += nucleotide_change[char]
				seq23wt = new_res_sequence[0:52].lower()
				mes23lowdiff, mes23highdiff = mes23highdiff - 1, mes23lowdiff + 1
				res_sequence = new_res_sequence
			else:
				seq23wt = res_sequence[0:52].lower()

			if "del" in hgvs:
				seq23mut = "{}{}".format(res_sequence[0:mes23lowdiff], res_sequence[mes23lowdiff+indel_length:len(res_sequence)]).lower()
			elif "ins" in hgvs:
				seq23wt += res_sequence[52:52+indel_length].lower()
				seq23mut = "{}{}{}".format(res_sequence[0:mes23lowdiff], indel_bp, res_sequence[mes23lowdiff:mes23lowdiff+mes23highdiff-len(indel_bp)]).lower()
			elif "dup" in hgvs:
				seq23wt += res_sequence[52:52+indel_length].lower()
				seq23mut = "{}{}{}".format(res_sequence[0:mes23lowdiff], indel_bp, res_sequence[mes23lowdiff:mes23lowdiff+mes23highdiff-len(indel_bp)]).lower()
			else:
				seq23mut = "{}{}{}".format(seq23wt[0:mes23lowdiff], mut, seq23wt[mes23lowdiff+1:mes23lowdiff+mes23highdiff]).lower()
			# print(res_sequence.lower(), mes23lowdiff, mes23highdiff, indel_length, indel_bp)
			# print(seq23wt, len(seq23wt))
			# print(seq23mut, len(seq23mut), "\n")
			return create_mes_dict(seq23wt, seq23mut, 3, strand)
		else:
			print("Something went wrong with 3' maxentscan. Please re-enter your coordinates.")
			return {}
	except Exception as e:
		print("Could not run 3' maxentscan: ", e)
		wt3mes = [False]
		mut3mes = [False]
		return {}

# determine aso amenability using metrics
def determine_amenability(mes_analyses, spliceai_analyses, hgvs, consequence, sift_score, polyphen_score, strand):

	sift_threshold = 0.05
	polyphen_threshold = 0.85
	deep_intron_threshold = 50 # Taken from the SpliceAI paper
	aso_amenability = ""
	intron, intron_depth = False, ""

	# Analyze the variant - in an intron?
	if "+" in hgvs or "-" in hgvs:
		intron = True
		if "+" in hgvs:
			for char in hgvs.split("+")[1]:
				if char.isnumeric():
					intron_depth += char
		elif "-" in hgvs:
			intron_depth += "-"
			for char in hgvs.split("-")[1]:
				if char.isnumeric():
					intron_depth += char	
		else:
			intron_depth = "0"
	else:
		intron_depth = "0"
	 # intron_depth = int(intron_depth) * strand
	intron_depth = int(intron_depth)

	# Step 1: First check for stop gain / frameshift
	if consequence == "stop_gained" or "stop" in consequence:
		aso_amenability = "It's unlikely this variant is ASO-amenable because it is a nonsense mutation."
		return aso_amenability
	elif consequence == "frameshift_variant" or "frameshift" in consequence:
		aso_amenability = "It's unlikely this variant is ASO-amenable because it is a frameshift mutation."
		return aso_amenability

	# # Step 2: Then first check pathogenicity
	# if type(sift_score) is not list and type(polyphen_score) is not list: 
	# 	if sift_score <= sift_threshold and polyphen_score >= polyphen_threshold:
	# 		aso_amenability = "It's very unlikely this variant is ASO-amenable because SIFT and PolyPhen scores predict that the mutation is deleterious."
	# 	# Exonic - benign. Check to see if there's a new site created nearby
	# 	elif sift_score >= sift_threshold and polyphen_score <= polyphen_threshold:
	# 		mes_count, combined_count, likely_count = 0, 0, 0	# found in mes, spliceai, and likely vs possible count
	# 		for entry in mes_analyses:
	# 			if int(mes_analyses[entry][2]) == 1:
	# 				mes_count += 1	# MES finds a site being strengthened
	# 				if mes_analyses[entry][0] == 2 or mes_analyses[entry][0] == 1 or float(entry.split(":")[3]) > 3.00:
	# 					likely_count += 1	# considered likely by MES
	# 			for spliceai_key in spliceai_analyses:
	# 				if float(entry.split(":")[0]) == (strand * float(spliceai_key.split(":")[0])):	# if splice ai finds an entry at the same position
	# 					if mes_analyses[entry][2] == 1 and spliceai_analyses[spliceai_key][1] == 1 and int(mes_analyses[entry][1]) == int(spliceai_analyses[spliceai_key][0]):	# both SpliceAI and MES indicate strengthening, and referring to same site
	# 						combined_count += 1
	# 		if mes_count == 0 and combined_count == 0:
	# 			aso_amenability = "It's very unlikely this variant is ASO-amenable. SIFT and PolyPhen scores predict that this mutation is benign/tolerated, but MES and SpliceAI analysis both do not predict any new sites being created/strengthened."
	# 		elif mes_count > 0 and combined_count == 0:
	# 			if likely_count == 0:
	# 				aso_amenability = "It's possible that this variant is ASO-amenable. SIFT and PolyPhen scores predict that this mutation is benign/tolerated, but only MES predicts a low probability that a new splice site is being created/strengthened."
	# 			else:
	# 				aso_amenability = "It's likely that this variant is ASO-amenable. SIFT and PolyPhen scores predict that this mutation is benign/tolerated, and MES alone predicts a moderate/high probability that a new splice site is being created/strengthened."
	# 		else:
	# 			aso_amenability = "It's very likely that this variant is ASO-amenable. SIFT and PolyPhen scores predict that this mutation is benign/tolerated, and MES and SpliceAI both predict that new sites are being created/strengthened."
	# 	else:	# SIFT/Polyphen conflict - return unsure
	# 		aso_amenability = "We're unsure if this variant is ASO-amenable. SIFT and PolyPhen scores have conflicting reports about the pathogenicity of the mutation."
	
	# Step 2: Then first check pathogenicity
	if type(sift_score) is not list and type(polyphen_score) is not list:
		#emsherr edit - don't bother to check SIFT/PolyPhen if synonymous
		if consequence != "synonymous_variant":
			if sift_score <= sift_threshold and polyphen_score >= polyphen_threshold:
				coding_impact = "SIFT and PolyPhen scores predict that the mutation is deleterious."
			# Exonic - benign. Check to see if there's a new site created nearby
			elif sift_score >= sift_threshold and polyphen_score <= polyphen_threshold:
				coding_impact = "SIFT and PolyPhen scores predict that this mutation is benign/tolerated."
			else:
				coding_impact = "SIFT and PolyPhen scores have conflicting reports about the pathogenicity of this variant."
		mes_count, combined_count, likely_count = 0, 0, 0	# found in mes, spliceai, and likely vs possible count
		#return mes_count
		for entry in mes_analyses:
			if int(mes_analyses[entry][2]) == 1:
				mes_count += 1	# MES finds a site being strengthened
				#return mes_count
				if mes_analyses[entry][0] == 2 or mes_analyses[entry][0] == 1 or float(entry.split(":")[3]) > 3.00:
					likely_count += 1	# considered likely by MES
			for spliceai_key in spliceai_analyses:
				if float(entry.split(":")[0]) == (strand * float(spliceai_key.split(":")[0])):	# if splice ai finds an entry at the same position
					if mes_analyses[entry][2] == 1 and spliceai_analyses[spliceai_key][1] == 1 and int(mes_analyses[entry][1]) == int(spliceai_analyses[spliceai_key][0]):	# both SpliceAI and MES indicate strengthening, and referring to same site
						combined_count += 1
						#return combined_count

		#emsherr edit
		if coding_impact == "SIFT and PolyPhen scores predict that the mutation is deleterious.":
			if mes_count == 0 and combined_count == 0:
				aso_amenability = "It's unlikely this variant is ASO-amenable. MES and SpliceAI analysis both do not predict any new sites being created/strengthened. {}".format(coding_impact)
			elif mes_count > 0 and combined_count == 0:
				if likely_count == 0:
					aso_amenability = "It's unlikely that this variant is ASO-amenable. Between MES and SpliceAI, only MES predicts a low probability that a new splice site is being created/strengthened and {}".format(coding_impact)
				else:
					aso_amenability = "It's unlikely that this variant is ASO-amenable. MES alone predicts a moderate/high probability that a new splice site is being created/strengthened and {}".format(coding_impact)
			else:
				aso_amenability = "It's unlikely that this variant is ASO-amenable. MES and SpliceAI both predict that new sites are being created/strengthened, but {}".format(coding_impact)

		else:
			if mes_count == 0 and combined_count == 0:
				aso_amenability = "It's unlikely this variant is ASO-amenable. MES and SpliceAI analysis both do not predict any new sites being created/strengthened. {}".format(coding_impact)
			elif mes_count > 0 and combined_count == 0:
				if likely_count == 0:
					aso_amenability = "It's possible that this variant is ASO-amenable. Between MES and SpliceAI, only MES predicts a low probability that a new splice site is being created/strengthened. {}".format(coding_impact)
				else:
					aso_amenability = "This variant is probably ASO-amenable. MES alone predicts a moderate/high probability that a new splice site is being created/strengthened. {}".format(coding_impact)
			else:
				aso_amenability = "This variant is probably ASO-amenable. MES and SpliceAI both predict that new sites are being created/strengthened. {}".format(coding_impact)


		#emsherr edit - check canonical here
		try:
			first_key = str(next(iter(mes_analyses)))
		except Exception as e:
			first_key = None
		if len(mes_analyses) == 1 and first_key and (len(first_key.split(":")) == 3 or len(first_key.split(":")) != 4):	# if VEP's MES results
			# Try to find canonical impact through SpliceAI results
			analyzed_canonical = False
			for spliceai_key in spliceai_analyses:
				if (strand * float(spliceai_key.split(":")[0])) == float(intron_depth):	# found canonical impact in spliceai results
					analyzed_canonical = True
					if spliceai_analyses[spliceai_key][1] == 0:	# if canonical site is found to be weakened, check MES result or other spliceai results to see if there is a site being strengthened nearby

						# first check if is destroyed - use Jinkuk's criteria to see if spliceAI returns >0.5 prob that it is weakened
						if float(spliceai_key.split(":")[1]) > 0.5:
							aso_amenability += " In addition, SpliceAI analysis predicts that the canonical splice site is most likely destroyed."
							break
						#emsherr add exon skipping?
						mes_strengthened, spliceai_strengthened = False, False
						if mes_analyses[first_key][2] == 1 and (mes_analyses[first_key][0] == 2 or mes_analyses[first_key][0] == 1):
							mes_strengthened = True
						for second_run in spliceai_analyses:
							if spliceai_analyses[second_run][1] == 1:	# another result found to be strengthened in spliceai
								spliceai_strengthened = True
						if mes_strengthened and spliceai_strengthened:
							aso_amenability += " In addition, SpliceAI analysis predicts that the canonical splice site is being weakened, while MES/SpliceAI results predict that a new splice site is being created/strengthened."
						elif not mes_strengthened and not spliceai_strengthened:
							#emsherrskipping
							aso_amenability += " In addition, SpliceAI analysis predicts that the canonical splice site is being weakened, but MES/SpliceAI results do not predict that a new splice site is being created/strengthened. It's possible that this mutation is weakening exon definition and may lead to exon skipping/aberrant splice patterns."
						else:
							aso_amenability += " In addition, SpliceAI analysis predicts that the canonical splice site is being weakened, but MES/SpliceAI do not agree that a new splice site is being created/strengthened."
					else:
						aso_amenability += " In addition, SpliceAI analysis predicts that the canonical splice site is not being weakened."

		
		# If we have MES results from Perl, check if the canonical splice site is still intact - look for a MES signal at the canonical splice site
		analyzed_canonical = False
		for entry in mes_analyses:
			pos = entry.split(":")[0]
			if abs(int(float(pos))) == abs(intron_depth) and intron_depth * -1 == int(pos):	# need to be opposite signs and equal position to assess the canonical splite site
				analyzed_canonical = True
				#emsherr update weakened not destroyed criteria
				#if mes_analyses[entry][2] == 0 or float(entry.split(":")[3]) < 4.00:	# The canonical site is predicted to be weakened
				if mes_analyses[entry][2] == 0 or abs(float(entry.split(":")[3])) < 0.7*float(entry.split(":")[1]):	# The canonical site is predicted to be weakened


					# destroyed - not amenable
					#emsherr
					#if float(entry.split(":")[2]) <= 3.00:
					if float(entry.split(":")[2]) < 2.00:
						aso_amenability += " Using MES, we analyzed that the canonical splice site is likely destroyed."
						break

					#emsherr try skipping here??
					mes_found = False
					beyond_results, within_results = [], [] 	# mapping of beyond/within 10: whether MES/SpliceAI support it
					splice_beyond, splice_within = [], []

					# Now check if there's a new site being created/strengthened
					for second_run in mes_analyses:	
						pos_gain = second_run.split(":")[0]
						if mes_analyses[second_run][2] == 1:
							direction = False
							if (int(pos_gain) + intron_depth >= 10):	# if a new site is being created/strengthened beyond 10 bp of the canonical splice site
								# CAN ALSO INPUT A PROBABILITY THRESHOLD HERE - IF MES SEES A HIGH/MODERATE PROBABILITY, ADD TO LIST
								beyond_results.append(entry)
							else:	# new site being created within 10 bp of the canonical splice site
								direction = True
								within_results.append(entry)

							# Look through SpliceAI to see if they confirm MES results
							counter = 0
							for spliceai_key in spliceai_analyses:
								# compare that the positions are the same and that they refer to the same splicing effect (strengthened/weakened vs increased/decreased)
								if float(spliceai_key.split(":")[0]) == float(pos) and int(spliceai_analyses[spliceai_key][1]) == int(mes_analyses[entry][2]):	# confirming canonical site
									counter += 1
								elif float(spliceai_key.split(":")[0]) == float(pos_gain) and int(spliceai_analyses[spliceai_key][1]) == int(mes_analyses[second_run][2]):	# confirming new splice site strengthening
									counter += 1
							if counter == 2:
								if direction:
									splice_within.append(spliceai_key)
								else:
									splice_beyond.append(spliceai_key)
							mes_found = True

					# Some new splice site was found - create the prediction
					if mes_found:
						if len(beyond_results) == 0 and len(within_results) > 0:
							aso_amenability = "It's unlikely this variant is ASO-amenable. Using MES, we analyzed that the canonical splice site is being weakened, while {} new splice site(s) are being created/strengthened within 10bp of the canonical site.".format(str(len(within_results)))
							if len(splice_within) > 0 and len(splice_within) != len(within_results):
								aso_amenability += " SpliceAI agrees with some of these assessments."
							elif len(splice_within) == len(within_results):
								aso_amenability += " SpliceAI agrees with all of those assessments."
						elif len(beyond_results) > 0 and len(within_results) == 0:
							aso_amenability = "This variant is probably ASO-amenable. Using MES, we analyzed that the canonical splice site is being weakened, while {} new splice site(s) are being created/strengthened beyond 10bp of the canonical site.".format(str(len(beyond_results)))
							if len(splice_beyond) > 0 and len(splice_beyond) != len(beyond_results):
								aso_amenability += " SpliceAI agrees with some of these assessments."		
							elif len(splice_beyond) == len(beyond_results):
								aso_amenability += " SpliceAI agrees with all of those assessments."
						elif len(beyond_results) > 0 and len(within_results) > 0:
							aso_amenability = "This variant is probably ASO-amenable. Using MES, we analyzed that the canonical splice site is being weakened, while {} new splice site(s) are being created/strengthened beyond 10bp of the canonical site and {} new splice site(s) are being created/strengthened within 10bp of the canonical site.".format(str(len(beyond_results)), str(len(within_results)))
							if len(splice_beyond) > 0 and len(splice_within) > 0 and (len(splice_beyond) != len(beyond_results) or len(splice_within) != len(within_results)):
								aso_amenability += " SpliceAI agrees with some of these assessments."		
							elif len(splice_beyond) == len(beyond_results) and len(splice_within) == len(within_results):
								aso_amenability += " SpliceAI agrees with all of those assessments."
						else:	# nothing found - aberrant splicing
							#emsherrskipping
							aso_amenability = "It's possible this variant is ASO-amenable. Using MES/SpliceAI, we analyzed that the canonical splice site is being weakened, while they do not predict a new splice site being created/strengthened nearby. It's possible that this mutation is weakening exon definition and may lead to exon skipping/aberrant splice patterns."

					# No new splice sites found using MES - try SpliceAI
					else:
						spliceai_found, spliceai_pos = False, ""
						for spliceai_key in spliceai_analyses:
							# if spliceai is referring to an increase
							if int(spliceai_analyses[spliceai_key][1]) == 1:
								spliceai_found, spliceai_pos = True, spliceai_key.split(":")[0]
						if not spliceai_found:
							#emsherrskipping
							aso_amenability = "It's possible this variant is ASO-amenable. Using MES/SpliceAI, we analyzed that the canonical splice site is being weakened, while they do not predict a new splice site being created/strengthened nearby. It's possible that this mutation is weakening exon definition and may lead to exon skipping/aberrant splice patterns."
						else:
							if float(spliceai_pos) >= -10 and float(spliceai_pos) <= 10:
								aso_amenability = "It's unlikely this variant is ASO-amenable. Using MES/SpliceAI, we analyzed that the canonical splice site is weakened, but only SpliceAI found that a new splice site was being created/strengthened within 10bp of the variant."
							else:
								aso_amenability = "This variant is probably ASO-amenable. Using MES/SpliceAI, we analyzed that the canonical splice site is weakened, but only SpliceAI found that a new splice site was being created/strengthened beyond 10bp of the variant."
				
					#what I had before
					#else:
						#aso_amenability += " In addition, MES predicts that the canonical splice site is being weakened."

				# canonical site is not weakened - check to see if it's destroyed 
				else:
					if float(entry.split(":")[2]) <= 4.00:
						aso_amenability += " Using MES, we analyzed that the canonical splice site is likely destroyed."
						break
					else:
						aso_amenability += " In addition, MES predicts that the canonical splice site is not being weakened."

	#emsherr edit - check synonymous variants without SIFT/PolyPhen	
	elif intron == False:
		mes_count, combined_count, likely_count = 0, 0, 0	# found in mes, spliceai, and likely vs possible count
		for entry in mes_analyses:
			if int(mes_analyses[entry][2]) == 1:
				mes_count += 1	# MES finds a site being strengthened
				#return mes_count
				if mes_analyses[entry][0] == 2 or mes_analyses[entry][0] == 1 or float(entry.split(":")[3]) > 3.00:
					likely_count += 1	# considered likely by MES
			for spliceai_key in spliceai_analyses:
				if float(entry.split(":")[0]) == (strand * float(spliceai_key.split(":")[0])):	# if splice ai finds an entry at the same position
					if mes_analyses[entry][2] == 1 and spliceai_analyses[spliceai_key][1] == 1 and int(mes_analyses[entry][1]) == int(spliceai_analyses[spliceai_key][0]):	# both SpliceAI and MES indicate strengthening, and referring to same site
						combined_count += 1
						#return combined_count

		#emsherr edit
		if mes_count == 0 and combined_count == 0:
			aso_amenability = "It's unlikely this variant is ASO-amenable. MES and SpliceAI analysis both do not predict any new sites being created/strengthened."
			#return aso_amenability
		elif mes_count > 0 and combined_count == 0:
			if likely_count == 0:
				aso_amenability = "It's possible that this variant is ASO-amenable. Between MES and SpliceAI, only MES predicts a low probability that a new splice site is being created/strengthened."
				#return aso_amenability
			else:
				aso_amenability = "This variant is probably ASO-amenable. MES alone predicts a moderate/high probability that a new splice site is being created/strengthened."
				#return aso_amenability
		else:
			aso_amenability = "This variant is probably ASO-amenable. MES and SpliceAI both predict that new sites are being created/strengthened."
		#return aso_amenability


		#emsherr edit - check canonical here
		try:
			first_key = str(next(iter(mes_analyses)))
		except Exception as e:
			first_key = None
		if len(mes_analyses) == 1 and first_key and (len(first_key.split(":")) == 3 or len(first_key.split(":")) != 4):	# if VEP's MES results
			# Try to find canonical impact through SpliceAI results
			analyzed_canonical = False
			for spliceai_key in spliceai_analyses:
				if (strand * float(spliceai_key.split(":")[0])) == float(intron_depth):	# found canonical impact in spliceai results
					analyzed_canonical = True
					if spliceai_analyses[spliceai_key][1] == 0:	# if canonical site is found to be weakened, check MES result or other spliceai results to see if there is a site being strengthened nearby

						# first check if is destroyed - use Jinkuk's criteria to see if spliceAI returns >0.5 prob that it is weakened
						if float(spliceai_key.split(":")[1]) > 0.5:
							aso_amenability += " In addition, SpliceAI analysis predicts that the canonical splice site is most likely destroyed."
							break
						#emsherr add exon skipping?
						mes_strengthened, spliceai_strengthened = False, False
						if mes_analyses[first_key][2] == 1 and (mes_analyses[first_key][0] == 2 or mes_analyses[first_key][0] == 1):
							mes_strengthened = True
						for second_run in spliceai_analyses:
							if spliceai_analyses[second_run][1] == 1:	# another result found to be strengthened in spliceai
								spliceai_strengthened = True
						if mes_strengthened and spliceai_strengthened:
							aso_amenability += " In addition, SpliceAI analysis predicts that the canonical splice site is being weakened, while MES/SpliceAI results predict that a new splice site is being created/strengthened."
						elif not mes_strengthened and not spliceai_strengthened:
							#emsherrskipping
							aso_amenability += " In addition, SpliceAI analysis predicts that the canonical splice site is being weakened, but MES/SpliceAI results do not predict that a new splice site is being created/strengthened. It's possible that this mutation is weakening exon definition and may lead to exon skipping/aberrant splice patterns."
						else:
							aso_amenability += " In addition, SpliceAI analysis predicts that the canonical splice site is being weakened, but MES/SpliceAI do not agree that a new splice site is being created/strengthened."
					else:
						aso_amenability += " In addition, SpliceAI analysis predicts that the canonical splice site is not being weakened."


		# If we have MES results from Perl, check if the canonical splice site is still intact - look for a MES signal at the canonical splice site
		analyzed_canonical = False
		for entry in mes_analyses:
			pos = entry.split(":")[0]
			if abs(int(float(pos))) == abs(intron_depth) and intron_depth * -1 == int(pos):	# need to be opposite signs and equal position to assess the canonical splite site
				analyzed_canonical = True
				#emsherr edit, changing criteria 
				#if mes_analyses[entry][2] == 0 or float(entry.split(":")[3]) < 4.00:	# The canonical site is predicted to be weakened
				if mes_analyses[entry][2] == 0 or abs(float(entry.split(":")[3])) < 0.7*float(entry.split(":")[1]):

					# destroyed - not amenable
					#emsherr
					#if float(entry.split(":")[2]) <= 3.00:
					if float(entry.split(":")[2]) < 2.00:
						aso_amenability += " Using MES, we analyzed that the canonical splice site is likely destroyed."
						break

					#emsherr try skipping here??
					mes_found = False
					beyond_results, within_results = [], [] 	# mapping of beyond/within 10: whether MES/SpliceAI support it
					splice_beyond, splice_within = [], []

					# Now check if there's a new site being created/strengthened
					for second_run in mes_analyses:	
						pos_gain = second_run.split(":")[0]
						if mes_analyses[second_run][2] == 1:
							direction = False
							if (int(pos_gain) + intron_depth >= 10):	# if a new site is being created/strengthened beyond 10 bp of the canonical splice site
								# CAN ALSO INPUT A PROBABILITY THRESHOLD HERE - IF MES SEES A HIGH/MODERATE PROBABILITY, ADD TO LIST
								beyond_results.append(entry)
							else:	# new site being created within 10 bp of the canonical splice site
								direction = True
								within_results.append(entry)

							# Look through SpliceAI to see if they confirm MES results
							counter = 0
							for spliceai_key in spliceai_analyses:
								# compare that the positions are the same and that they refer to the same splicing effect (strengthened/weakened vs increased/decreased)
								if float(spliceai_key.split(":")[0]) == float(pos) and int(spliceai_analyses[spliceai_key][1]) == int(mes_analyses[entry][2]):	# confirming canonical site
									counter += 1
								elif float(spliceai_key.split(":")[0]) == float(pos_gain) and int(spliceai_analyses[spliceai_key][1]) == int(mes_analyses[second_run][2]):	# confirming new splice site strengthening
									counter += 1
							if counter == 2:
								if direction:
									splice_within.append(spliceai_key)
								else:
									splice_beyond.append(spliceai_key)
							mes_found = True

					# Some new splice site was found - create the prediction
					if mes_found:
						if len(beyond_results) == 0 and len(within_results) > 0:
							aso_amenability = "It's unlikely this variant is ASO-amenable. Using MES, we analyzed that the canonical splice site is being weakened, while {} new splice site(s) are being created/strengthened within 10bp of the canonical site.".format(str(len(within_results)))
							if len(splice_within) > 0 and len(splice_within) != len(within_results):
								aso_amenability += " SpliceAI agrees with some of these assessments."
							elif len(splice_within) == len(within_results):
								aso_amenability += " SpliceAI agrees with all of those assessments."
						elif len(beyond_results) > 0 and len(within_results) == 0:
							aso_amenability = "This variant is probably ASO-amenable. Using MES, we analyzed that the canonical splice site is being weakened, while {} new splice site(s) are being created/strengthened beyond 10bp of the canonical site.".format(str(len(beyond_results)))
							if len(splice_beyond) > 0 and len(splice_beyond) != len(beyond_results):
								aso_amenability += " SpliceAI agrees with some of these assessments."		
							elif len(splice_beyond) == len(beyond_results):
								aso_amenability += " SpliceAI agrees with all of those assessments."
						elif len(beyond_results) > 0 and len(within_results) > 0:
							aso_amenability = "This variant is probably ASO-amenable. Using MES, we analyzed that the canonical splice site is being weakened, while {} new splice site(s) are being created/strengthened beyond 10bp of the canonical site and {} new splice site(s) are being created/strengthened within 10bp of the canonical site.".format(str(len(beyond_results)), str(len(within_results)))
							if len(splice_beyond) > 0 and len(splice_within) > 0 and (len(splice_beyond) != len(beyond_results) or len(splice_within) != len(within_results)):
								aso_amenability += " SpliceAI agrees with some of these assessments."		
							elif len(splice_beyond) == len(beyond_results) and len(splice_within) == len(within_results):
								aso_amenability += " SpliceAI agrees with all of those assessments."
						else:	# nothing found - aberrant splicing
							#emsherrskipping
							aso_amenability = "It's possible this variant is ASO-amenable. Using MES/SpliceAI, we analyzed that the canonical splice site is being weakened, while they do not predict a new splice site being created/strengthened nearby. It's possible that this mutation is weakening exon definition and may lead to exon skipping/aberrant splice patterns."

					# No new splice sites found using MES - try SpliceAI
					else:
						spliceai_found, spliceai_pos = False, ""
						for spliceai_key in spliceai_analyses:
							# if spliceai is referring to an increase
							if int(spliceai_analyses[spliceai_key][1]) == 1:
								spliceai_found, spliceai_pos = True, spliceai_key.split(":")[0]
						if not spliceai_found:
							#emsherrskipping
							aso_amenability = "It's possible this variant is ASO-amenable. Using MES/SpliceAI, we analyzed that the canonical splice site is being weakened, while they do not predict a new splice site being created/strengthened nearby. It's possible that this mutation is weakening exon definition and may lead to exon skipping/aberrant splice patterns."
						else:
							if float(spliceai_pos) >= -10 and float(spliceai_pos) <= 10:
								aso_amenability = "It's unlikely this variant is ASO-amenable. Using MES/SpliceAI, we analyzed that the canonical splice site is weakened, but only SpliceAI found that a new splice site was being created/strengthened within 10bp of the variant."
							else:
								aso_amenability = "This variant is probably ASO-amenable. Using MES/SpliceAI, we analyzed that the canonical splice site is weakened, but only SpliceAI found that a new splice site was being created/strengthened beyond 10bp of the variant."
				
					#what I had before
					#else:
						#aso_amenability += " In addition, MES predicts that the canonical splice site is being weakened."

				# canonical site is not weakened - check to see if it's destroyed 
				else:
					if float(entry.split(":")[2]) <= 4.00:
						aso_amenability += " Using MES, we analyzed that the canonical splice site is likely destroyed."
						break
					else:
						aso_amenability += " In addition, MES predicts that the canonical splice site is not being weakened."

	# Step 3: Go through MES/Splice AI for analysis
	else:
		if intron: # Intronic

			# First check if mes_analyses contains only VEP's MES results 
			try:
				first_key = str(next(iter(mes_analyses)))
			except Exception as e:
				first_key = None
			if len(mes_analyses) == 1 and first_key and (len(first_key.split(":")) == 3 or len(first_key.split(":")) != 4):	# if VEP's MES results
				# Try to find canonical impact through SpliceAI results
				analyzed_canonical = False
				for spliceai_key in spliceai_analyses:
					if (strand * float(spliceai_key.split(":")[0])) == float(intron_depth):	# found canonical impact in spliceai results
						analyzed_canonical = True
						if spliceai_analyses[spliceai_key][1] == 0:	# if canonical site is found to be weakened, check MES result or other spliceai results to see if there is a site being strengthened nearby

							# first check if is destroyed - use Jinkuk's criteria to see if spliceAI returns >0.5 prob that it is weakened
							if float(spliceai_key.split(":")[1]) > 0.5:
								aso_amenability = "It's unlikely that this variant is ASO-amenable. SpliceAI analysis predicts that the canonical splice site is most likely destroyed."
								break

							mes_strengthened, spliceai_strengthened = False, False
							if mes_analyses[first_key][2] == 1 and (mes_analyses[first_key][0] == 2 or mes_analyses[first_key][0] == 1):
								mes_strengthened = True
							for second_run in spliceai_analyses:
								if spliceai_analyses[second_run][1] == 1:	# another result found to be strengthened in spliceai
									spliceai_strengthened = True
							if mes_strengthened and spliceai_strengthened:
								aso_amenability = "It's possible that this variant is ASO-amenable. SpliceAI analysis predicts that the canonical splice site is being weakened, while MES/SpliceAI results predict that a new splice site is being created/strengthened."
							elif not mes_strengthened and not spliceai_strengthened:
								#emsherrskipping
								aso_amenability = "It's possible that this variant is ASO-amenable. SpliceAI analysis predicts that the canonical splice site is being weakened, but MES/SpliceAI results do not predict that a new splice site is being created/strengthened. It's possible that this mutation is weakening exon definition and may lead to exon skipping/aberrant splice patterns."
							else:
								aso_amenability = "It's unlikely that this variant is ASO-amenable. SpliceAI analysis predicts that the canonical splice site is being weakened, but MES/SpliceAI do not agree that a new splice site is being created/strengthened."
						else:
							aso_amenability = "It's unlikely that this variant is ASO-amenable. SpliceAI analysis predicts that the canonical splice site is not being weakened, and MES results are not extensive enough to make a prediction."
				
				# noting was found about the canonical site from SpliceAI
				if not analyzed_canonical:
					if abs(intron_depth) <= 10:	# Less than 10bp away from exon start/end, so unlikely
						#emsherr edit - check mes even if less than 10bp away
						if mes_analyses[first_key][2] == 1: # if found to be strengthened
							if mes_analyses[first_key][0] == 2 or mes_analyses[first_key][0] == 1 or float(first_key.split(":")[2]) < -3.00:	# high/moderate probability or delta < -3.00 (VEP uses a different diff calculation)
								aso_amenability = "It's possible that this variant is ASO-amenable, as VEP's MES results predict a high/moderate probability that a new site is being created/strengthened, but this variant is within 10bp of a canonical splice site."
							else:
								aso_amenability = "It's possible that this variant is ASO-amenable, as VEP's MES results predict a low probability that a new site is being created/strengthened, but this variant is within 10bp of a canonical splice site."
							for spliceai_key in spliceai_analyses:
								if spliceai_analyses[spliceai_key][1] == 1:	# if spliceAI also finds a site increasing nearby
									aso_amenability += " SpliceAI also predicts a site being strengthened nearby."
						else:
							aso_amenability = "It's unlikely that this variant is ASO-amenable, as this variant is within 10bp of a canonical splice site and no extensive MES results were found."
					else:
						if mes_analyses[first_key][2] == 1: # if found to be strengthened
							if mes_analyses[first_key][0] == 2 or mes_analyses[first_key][0] == 1 or float(first_key.split(":")[2]) < -3.00:	# high/moderate probability or delta < -3.00 (VEP uses a different diff calculation)
								aso_amenability = "This variant is probably ASO-amenable, as this variant is beyond 10bp of a canonical splice site and VEP's MES results predict a high/moderate probability that a new site is being created/strengthened."
							else:
								aso_amenability = "It's possible that this variant is ASO-amenable, as this variant is beyond 10bp of a canonical splice site, but VEP's MES results predict a low probability that a new site is being created/strengthened."
							for spliceai_key in spliceai_analyses:
								if spliceai_analyses[spliceai_key][1] == 1:	# if spliceAI also finds a site increasing nearby
									aso_amenability += " SpliceAI also predicts a site being strengthened nearby."
								break

				# Nothing found - refer back to user just in case
				if len(aso_amenability) == 0: 
					aso_amenability = "We're unsure if this variant is ASO-amenable. Please review MES/SpliceAI results for further information."
				return aso_amenability

			# If we have MES results from Perl, check if the canonical splice site is still intact - look for a MES signal at the canonical splice site
			analyzed_canonical = False
			for entry in mes_analyses:
				pos = entry.split(":")[0]
				if abs(int(float(pos))) == abs(intron_depth) and intron_depth * -1 == int(pos):	# need to be opposite signs and equal position to assess the canonical splite site
					analyzed_canonical = True
					#emsherr edit to change criteria
					#if mes_analyses[entry][2] == 0 or float(entry.split(":")[3]) < 4.00:	# The canonical site is predicted to be weakened
					if mes_analyses[entry][2] == 0 or abs(float(entry.split(":")[3])) < 0.7*float(entry.split(":")[1]):

						# destroyed - not amenable
						#emsherr edit
						#if float(entry.split(":")[2]) <= 3.00:
						if float(entry.split(":")[2]) < 2.00:
							aso_amenability = "It's unlikely this variant is ASO-amenable. Using MES, we analyzed that the canonical splice site is likely destroyed."
							break

						mes_found = False
						beyond_results, within_results = [], [] 	# mapping of beyond/within 10: whether MES/SpliceAI support it
						splice_beyond, splice_within = [], []

						# Now check if there's a new site being created/strengthened
						for second_run in mes_analyses:	
							pos_gain = second_run.split(":")[0]
							if mes_analyses[second_run][2] == 1:
								direction = False
								if (int(pos_gain) + intron_depth >= 10):	# if a new site is being created/strengthened beyond 10 bp of the canonical splice site
									# CAN ALSO INPUT A PROBABILITY THRESHOLD HERE - IF MES SEES A HIGH/MODERATE PROBABILITY, ADD TO LIST
									beyond_results.append(entry)
								else:	# new site being created within 10 bp of the canonical splice site
									direction = True
									within_results.append(entry)

								# Look through SpliceAI to see if they confirm MES results
								counter = 0
								for spliceai_key in spliceai_analyses:
									# compare that the positions are the same and that they refer to the same splicing effect (strengthened/weakened vs increased/decreased)
									if float(spliceai_key.split(":")[0]) == float(pos) and int(spliceai_analyses[spliceai_key][1]) == int(mes_analyses[entry][2]):	# confirming canonical site
										counter += 1
									elif float(spliceai_key.split(":")[0]) == float(pos_gain) and int(spliceai_analyses[spliceai_key][1]) == int(mes_analyses[second_run][2]):	# confirming new splice site strengthening
										counter += 1
								if counter == 2:
									if direction:
										splice_within.append(spliceai_key)
									else:
										splice_beyond.append(spliceai_key)
								mes_found = True

						# Some new splice site was found - create the prediction
						if mes_found:
							if len(beyond_results) == 0 and len(within_results) > 0:
								aso_amenability = "It's unlikely this variant is ASO-amenable. Using MES, we analyzed that the canonical splice site is being weakened, while {} new splice site(s) are being created/strengthened within 10bp of the canonical site.".format(str(len(within_results)))
								if len(splice_within) > 0 and len(splice_within) != len(within_results):
									aso_amenability += " SpliceAI agrees with some of these assessments."
								elif len(splice_within) == len(within_results):
									aso_amenability += " SpliceAI agrees with all of those assessments."
							elif len(beyond_results) > 0 and len(within_results) == 0:
								aso_amenability = "This variant is probably ASO-amenable. Using MES, we analyzed that the canonical splice site is being weakened, while {} new splice site(s) are being created/strengthened beyond 10bp of the canonical site.".format(str(len(beyond_results)))
								if len(splice_beyond) > 0 and len(splice_beyond) != len(beyond_results):
									aso_amenability += " SpliceAI agrees with some of these assessments."		
								elif len(splice_beyond) == len(beyond_results):
									aso_amenability += " SpliceAI agrees with all of those assessments."
							elif len(beyond_results) > 0 and len(within_results) > 0:
								aso_amenability = "This variant is probably ASO-amenable. Using MES, we analyzed that the canonical splice site is being weakened, while {} new splice site(s) are being created/strengthened beyond 10bp of the canonical site and {} new splice site(s) are being created/strengthened within 10bp of the canonical site.".format(str(len(beyond_results)), str(len(within_results)))
								if len(splice_beyond) > 0 and len(splice_within) > 0 and (len(splice_beyond) != len(beyond_results) or len(splice_within) != len(within_results)):
									aso_amenability += " SpliceAI agrees with some of these assessments."		
								elif len(splice_beyond) == len(beyond_results) and len(splice_within) == len(within_results):
									aso_amenability += " SpliceAI agrees with all of those assessments."
							else:	# nothing found - aberrant splicing
								#emsherrskipping
								aso_amenability = "It's possible this variant is ASO-amenable. Using MES/SpliceAI, we analyzed that the canonical splice site is being weakened, while they do not predict a new splice site being created/strengthened nearby. It's possible that this mutation is weakening exon definition and may lead to exon skipping/aberrant splice patterns."

						# No new splice sites found using MES - try SpliceAI
						else:
							spliceai_found, spliceai_pos = False, ""
							for spliceai_key in spliceai_analyses:
								# if spliceai is referring to an increase
								if int(spliceai_analyses[spliceai_key][1]) == 1:
									spliceai_found, spliceai_pos = True, spliceai_key.split(":")[0]
							if not spliceai_found:
								#emsherrskipping
								aso_amenability = "It's possible this variant is ASO-amenable. Using MES/SpliceAI, we analyzed that the canonical splice site is being weakened, while they do not predict a new splice site being created/strengthened nearby. It's possible that this mutation is weakening exon definition and may lead to exon skipping/aberrant splice patterns."
							else:
								if float(spliceai_pos) >= -10 and float(spliceai_pos) <= 10:
									aso_amenability = "It's unlikely this variant is ASO-amenable. Using MES/SpliceAI, we analyzed that the canonical splice site is weakened, but only SpliceAI found that a new splice site was being created/strengthened within 10bp of the variant."
								else:
									aso_amenability = "This variant is probably ASO-amenable. Using MES/SpliceAI, we analyzed that the canonical splice site is weakened, but only SpliceAI found that a new splice site was being created/strengthened beyond 10bp of the variant."
					
					# canonical site is not weakened - check to see if it's destroyed 
					else:
						if float(entry.split(":")[2]) <= 4.00:
							aso_amenability = "It's unlikely this variant is ASO-amenable. Using MES, we analyzed that the canonical splice site is likely destroyed."
							break
			
			# went through MES but no results on canonical splite site. See if it is deep intronic? What is the deep intronic threshold?
			if not analyzed_canonical and len(aso_amenability) == 0:
				if abs(intron_depth) >= deep_intron_threshold:	# deep intronic
					spliceai_agree = False
					for entry in mes_analyses:
						pos = entry.split(":")[0]
						if mes_analyses[entry][2] == 1:	# new site created/strengthened
							if mes_analyses[entry][0] == 2 or mes_analyses[entry][0] == 1:	# high/moderate prob
								aso_amenability = "This variant is probably ASO-amenable. Using MES, we analyzed that there is a splice site being created/strengthened nearby."
								for spliceai_key in spliceai_analyses:
									if float(spliceai_key.split(":")[0]) == float(pos) and spliceai_analyses[spliceai_key][0] == mes_analyses[entry][1] and spliceai_analyses[spliceai_key][1] == mes_analyses[entry][2]:
										aso_amenability += " SpliceAI agrees with this assessment."
								break
							else:	# low prob
								aso_amenability = "It's possible that this variant is ASO-amenable. Using MES, we analyzed that there is a low probability that a splice site is being created/strengthened nearby."
								for spliceai_key in spliceai_analyses:
									if spliceai_key.split(":")[0] == pos and spliceai_analyses[spliceai_key][0] == mes_analyses[entry][1] and spliceai_analyses[spliceai_key][1] == mes_analyses[entry][2]:
										aso_amenability += " SpliceAI agrees with this assessment."
					if len(aso_amenability) == 0:	# deep intronic but no site found to be created/strengthened
						aso_amenability = "It's possible that this variant is ASO-amenable. Using MES/SpliceAI, we did not find a new splice site being created/strengthened nearby. It's possible the variant is benign."
				else:	# not deep intronic - check spliceai for an impact on canonical site?
					canonical_found = False
					for spliceai_key in spliceai_analyses:
						pos = spliceai_key.split(":")[0]
						if float(pos) == float(abs(intron_depth)) and float(intron_depth * -1) == float(pos):	# found canonical splice site
							canonical_found = True
							if spliceai_analyses[spliceai_key][1] == 1:	# found that canonical splice site is likely being weakened
								for second_run in spliceai_analyses:	# look for any other results that indicate a strengthening/creation of a new site
									pos_gain = second_run.split(":")[0]
									if splice_analyses[second_run][1] == 1: # found another spliceai result that indicates a new site is being used
										if int(pos_gain) + intron_depth >= 10:
											aso_amenability = "This variant is probably ASO-amenable. Using SpliceAI, we analyzed that the canonical splite site is being weakened while a new splice site is being created/strengthened more than 10bp away from the original."
										else:
											aso_amenability = "It's unlikely this variant is ASO-amenable. Using SpliceAI, we analyzed that the canonical splite site is being weakened while a new splice site is being created/strengthened within 10bp of the original."
											break
					if not canonical_found:	# maybe between 10-50bp away?
						created_site_found = False
						for spliceai_key in spliceai_analyses:
							if spliceai_analyses[spliceai_key][1] == 1:	# found a site possibly being strengthened/created
								if abs(intron_depth) >= 10:		# beyond 10bp away of exon start/end
									# assess the probability given by spliceai - these are arbitrary thresholds - need to discuss
									prob, created_site_found = float(spliceai_key.split(":")[1]), True
									if prob < 0.5:
										aso_amenability = "It's possible this variant is ASO-amenable. Using SpliceAI, we analyzed that there is <50{} probability that a splice site is being created/strengthened nearby, but no information was found about the canonical site.".format("%")
									elif prob < 0.75:
										aso_amenability = "This variant is probably ASO-amenable. Using SpliceAI, we analyzed that there is 50-75{} probability that a splice site is being created/strengthened nearby, but no information was found about the canonical site.".format("%")
									else:
										aso_amenability = "It's possible this variant is ASO-amenable. Using SpliceAI, we analyzed that there is 75-100{} probability that a splice site being created/strengthened nearby, but no information was found about the canonical site.".format("%")
						if not created_site_found:
							aso_amenability = "It's unlikely this variant is ASO-amenable. Using SpliceAI and MES, we could not either determine an impact on the canonical splice site or find a splice site being created/strengthened nearby. Please review the MES and SpliceAI results for further information."
	
	# still 0 - no usable results from SpliceAI or MES, refer back to user
	# maybe evaluate the intron depth?
	if len(aso_amenability) == 0: 
		if intron and abs(intron_depth) < 10:
			aso_amenability = "It's unlikely this variant is ASO-amenable, as it is extremely close to the original splice site and no other information was found. Please review MES/SpliceAI results for further information."
		else:	# no clue
			aso_amenability = "We're unsure if this variant is ASO-amenable. Please review MES/SpliceAI results for further information."

	return aso_amenability
	#return [aso_amenability, mes_count, combined_count]

# determine aso amenabiltiy using Jinkuk's thresholds
def determine_amenability_jinkuk(mes_analyses, spliceai_analyses, hgvs, consequence, sift_score, polyphen_score, strand):
	sift_threshold = 0.05
	polyphen_threshold = 0.85
	deep_intron_threshold = 50	# Taken from the SpliceAI paper
	aso_amenability = ""
	intron, intron_depth = False, ""
	coding_impact_numerical, splicing_impact_numerical = 0, 3

	# Analyze the variant - in an intron?
	if "+" in hgvs or "-" in hgvs:
		intron = True
		if "+" in hgvs:
			for char in hgvs.split("+")[1]:
				if char.isnumeric():
					intron_depth += char
		elif "-" in hgvs:
			intron_depth += "-"
			for char in hgvs.split("-")[1]:
				if char.isnumeric():
					intron_depth += char	
		else:
			intron_depth = "0"
	else:
		intron_depth = "0"
	 # intron_depth = int(intron_depth) * strand
	intron_depth = int(intron_depth)

	# Step 1: First check for stop gain / frameshift
	if consequence == "stop_gained" or "stop" in consequence:
		aso_amenability = "It's very unlikely this variant is ASO-amenable because it is a nonsense mutation."
		return aso_amenability
	elif consequence == "frameshift_variant" or "frameshift" in consequence:
		aso_amenability = "It's very unlikely this variant is ASO-amenable because it is a frameshift mutation."
		return aso_amenability

	# Step 2: Analyze coding damage
	if type(sift_score) is not list and type(polyphen_score) is not list: 
		if sift_score <= sift_threshold and polyphen_score >= polyphen_threshold:	# complete damage
			coding_impact_numerical = 1
		elif sift_score >= sift_threshold and polyphen_score <= polyphen_threshold:	# no or little
			coding_impact_numerical = 3
		else:	# all else
			coding_impact_numerical = 2

	# Step 3: Analyze splicing damage
	for entry in spliceai_analyses:
		if len(spliceai_analyses) == 1 and spliceai_analyses[entry][1] == 0 and float(entry.split(":")[1]) < 0.02:	# no or little
			new_impact = 3
			if new_impact <= splicing_impact_numerical:
				splicing_impact_numerical = new_impact
		elif spliceai_analyses[entry][1] == 0 and float(entry.split(":")[1]) >= 0.02 and float(entry.split(":")[1]) < 0.1:	# moderate
			for m_entry in mes_analyses:
				if float(m_entry.split(":")[1]) >= 2 and float(m_entry.split(":")[1]) >= (0.03 * float(m_entry.split(":")[0])):
					new_impact = 2
					if new_impact <= slicing_impact_numerical:
						splicing_impact_numerical = new_impact
		elif spliceai_analyses[entry][1] == 0 and float(entry.split(":")[1]) >= 0.1:	# strong/incomplete and complete damage
			for m_entry in mes_analyses:
				if float(m_entry.split(":")[1]) >= 2 and float(m_entry.split(":")[1]) >= (0.03 * float(m_entry.split(":")[0])):	# strong/incomplete
					new_impact = 1
					if new_impact <= slicing_impact_numerical:
						splicing_impact_numerical = new_impact
				elif float(m_entry.split(":")[1]) < 2 or float(m_entry.split(":")[1]) < (0.03 * float(m_entry.split(":")[0])):	# complete damage
					new_impact = 1
					if new_impact <= slicing_impact_numerical:
						splicing_impact_numerical = new_impact

# refseqID: NM_003159
# ensemblID: ENST00000379996 for cdkl5
def get_output(hgvs, wt="", mut="", transcript=""):
	"""Call tools to get output"""
	# Run VEP on the variant

	try:
		server = "https://rest.ensembl.org"
		get = "/vep/human/hgvs/" + hgvs
		params = {"SpliceAI": True, "MaxEntScan": True, "CADD": True}
		if transcript:
			params["transcript_id"] = transcript
		r = requests.get(server + get, params=params, headers={"Content-Type": "application/json"}, verify=False)

		if not r.ok:
			r.raise_for_status()
			return None

		decoded = r.json()
		res = decoded[0]
		print(res)
	except Exception as e:
		print("Could not run VEP: ", e)
		return None

	# adjust for strand
	strand, indel_length, indel_bp = 1, 0, ""
	for key in res:
		# print(key, res[key])
		if (key == "strand" or "strand" in key) and (res[key] == -1):
			strand = -1
		if "del" in hgvs and (key == "allele_string" or "allele_string" in key):
			indel_info = res[key].strip().split("/")
			indel_length, indel_bp = len(indel_info[0]), str(indel_info[0])
		elif "ins" in hgvs and (key == "allele_string" or "allele_string" in key):
			indel_info = res[key].strip().split("/")
			indel_length, indel_bp = len(indel_info[1]), str(indel_info[1])
		elif "dup" in hgvs and (key == "allele_string" or "allele_string" in key):
			indel_info = res[key].strip().split("/")
			indel_length, indel_bp = len(indel_info[1]), str(indel_info[1])
	if strand == -1 and len(wt) != 0 and len(mut) != 0:
		splice_wt, splice_mut = nucleotide_change[wt], nucleotide_change[mut]
	else:
		splice_wt, splice_mut = wt, mut

	# find SIFT and Polyphen scores
	try:
		sift_score = [False]
		for item in res["transcript_consequences"]:
			if "sift_score" in item:
				sift_score = item["sift_score"]
				break
	except Exception as e:
		print("Could not retrieve sift_score: ", e)
		sift_score = [False]
	try:
		polyphen_score = [False]
		for item in res["transcript_consequences"]:
			if "polyphen_score" in item:
				polyphen_score = item["polyphen_score"]
				break
	except Exception as e:
		print("Could not retrieve polyphen_score: ", e)
		polyphen_score = [False]

	# find CADD score
	try:
		cadd_score = [False]
		for item in res["transcript_consequences"]:
			if "cadd_phred" in item:
				cadd_score = item["cadd_phred"]
				break
	except Exception as e:
		print("Could not retrieve cadd_score: ", e)
		cadd_score = [False]

	# find these scores (maxentscan, spliceAI)
	try:
		maxentscan_ref, ref_found = [False], False
		maxentscan_alt, alt_found = [False], False
		maxentscan_diff, diff_found = [False], False
		spliceai_pred, pred_found = [False], False
		for item in res["transcript_consequences"]:
			if "maxentscan_ref" in item and not ref_found:
				maxentscan_ref = item["maxentscan_ref"]
				ref_found = True
			if "maxentscan_alt" in item and not alt_found:
				maxentscan_alt = item["maxentscan_alt"]
				alt_found = True
			if "maxentscan_diff" in item and not diff_found:
				maxentscan_diff = item["maxentscan_diff"]
				diff_found = True
			#emsherr note - this is searching for wrong item, should be just "splicai" - ok because searching directly later 
			if "spliceai_pred" in item and not pred_found:
				spliceai_pred = item["spliceai_pred"]
				pred_found = True
	except Exception as e:
		print("Could not retrieve some MES score: ", e)

	# look for start and chr of mutation
	try:
		gen_start = res["start"]
		gen_start = str(int(gen_start) - 1)		# Adjust for maxentscan perl script
		if "seq_region_name" in res:
			gen_chr = res["seq_region_name"]
		elif "colocated_variants" in res:
			gen_chr = res["colocated_variants"][0]["seq_region_name"]
		else:
			gen_chr = "un"
	except Exception as e:
		print("Could not retrieve some needed location information: ", e)

	# look for most severe consequence
	try:
		consequence = "Not found"
		if "most_severe_consequence" in res:
			consequence = res["most_severe_consequence"]
	except Exception as e:
		print("Could not retrieve most severe consequence: ", e)

	# look for amino acid change
	try:
		aa_change, codons, cdna = "Not found", "not found", "Not found"
		for item in res["transcript_consequences"]:
			if "amino_acids" in item:
				aa_change = item["amino_acids"]
			if "codons" in item:
				codons = item["codons"]
			if "cdna_start" in item:
				cdna = item["cdna_start"]
	except Exception as e:
		print("Could not retrieve amino acid change: ", e)

	# create gnomad variant string 
	gnomad_variant = ""
	if gen_chr and gen_start and len(splice_wt) != 0 and len(splice_mut) != 0:
		gnomad_variant = "{}-{}-{}-{}".format(gen_chr, str(int(gen_start) + 1), splice_wt, splice_mut)

	# # get spliceai prediction from spliceai-lookup and mes scores
	spliceai_pred = spliceai(gnomad_variant)
	mes_dict_5 = mes5_runner(gen_start, gen_chr, mut, strand, hgvs, indel_length, indel_bp)
	mes_dict_3 = mes3_runner(gen_start, gen_chr, mut, strand, hgvs, indel_length, indel_bp)

	# Determine if this mutation is amenable
	coding_impact, splicing_impact, consequence, mes_analyses, spliceai_analyses = determine_impacts(gen_chr, gen_start, hgvs, consequence, sift_score, polyphen_score, maxentscan_ref, maxentscan_alt, maxentscan_diff, spliceai_pred, mes_dict_5, mes_dict_3, strand)
	
	# IMPT:
	# mes_analyses is a dictionary that maps key: "key,ref,alt,delta" to a tuple value: (x, y, z) where x=high/moderate/low prob, y=donor/acceptor, z=weakened/strengthened
	# spliceai_analyses is a dictionary that maps key: "pos_delta,prob_delta" to a tuple value: (x, y), where x=acceptor/donor, y=loss/gain
	print(mes_analyses, spliceai_analyses)
	for entry in mes_analyses:
		print(mes_analyses[entry][2])
	#print(determine_amenability(mes_analyses, spliceai_analyses, hgvs, consequence, sift_score, polyphen_score, strand))
	#print(determine_amenability(mes_analyses, spliceai_analyses, hgvs, consequence, sift_score, polyphen_score, strand)[1])
	#print(determine_amenability(mes_analyses, spliceai_analyses, hgvs, consequence, sift_score, polyphen_score, strand)[2])
	aso_prediction = determine_amenability(mes_analyses, spliceai_analyses, hgvs, consequence, sift_score, polyphen_score, strand)

	d = {
		"aso_prediction": aso_prediction,
		"coding_impact": coding_impact,
		"splicing_impact": splicing_impact,
		"sift_score": sift_score,
		"polyphen_score": polyphen_score,
		"maxentscan_ref": maxentscan_ref,
		"maxentscan_alt": maxentscan_alt,
		"maxentscan_diff": maxentscan_diff,
		"spliceai_pred": spliceai_pred,
		"gnomad_variant": gnomad_variant,
		"cadd_score": cadd_score,
		"strand": str(strand),
		"aa_change": aa_change,
		"codons": codons,
		"cdna": cdna
	}
	return d

def get_literature(hgvs):
	literature = []
	try:
		# Retrieve PubMed IDs
		server = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
		get = "esearch.fcgi?"
		params = {"db": "pubmed", "term":hgvs, "retmode":"json"}
		r = requests.get(server + get, params=params, headers={"Content-Type": "application/json"})

		if not r.ok:
			r.raise_for_status()
			return []
		res = r.json()
		id_list = res["esearchresult"]["idlist"]

		# go through ID list and get the title info
		for pubmed_id in id_list:
			try:
				# Retrieve titles
				server = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
				get = "esummary.fcgi?"
				params = {"db": "pubmed", "id":pubmed_id, "retmode":"json"}
				r = requests.get(server + get, params=params, headers={"Content-Type": "application/json"})

				if not r.ok:
					r.raise_for_status()
					return []
				res = r.json()
				literature.append((str(pubmed_id), str(res["result"][pubmed_id]["title"])))		
			except Exception as e:
				return []

	except Exception as e:
		return []

	return literature

def get_clinvar(hgvs):
	try:
		# Retrieve Clinvar ID
		server = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
		get = "esearch.fcgi?"
		params = {"db": "clinvar", "term":hgvs, "retmode":"json"}
		r = requests.get(server + get, params=params, headers={"Content-Type": "application/json"})
		if not r.ok:
			r.raise_for_status()
			return ""
		res = r.json()
		if int(res["esearchresult"]["count"]) >= 1:
			return str(res["esearchresult"]["idlist"][0])
		return ""
	except Exception as e:
		return ""


# refseqID: NM_003159
# ensemblID: ENST00000379996 for cdkl5
# get output given list of hgvs, wts, muts, transcript (?)
def get_output_list(hgvs, wt="", mut="", transcript=""):
	# clear carriage return
	hgvs = [i.strip() for i in hgvs]
	if wt == "":
		wt = [""] * len(hgvs)
	if mut == "":
		mut = [""] * len(hgvs)
	if transcript == "":
		transcript = [""] * len(hgvs)
	"""Call tools to get output"""
	# Run VEP on the variant
	# https://rest.ensembl.org/ POST
	try:
		server = "http://rest.ensembl.org"
		ext = "/vep/human/hgvs/"
		params = {"SpliceAI": True, "MaxEntScan": True, "CADD": True}
		if transcript:
			params["transcript_id"] = transcript
		headers = {"Content-Type": "application/json", "Accept": "application/json"}
		# format hgvs list for post request
		postData = json.dumps(hgvs)
		postData = '{ "hgvs_notations" : ' + postData + ' }'
		# print(postData)
		r = requests.post(server+ext, params=params, headers=headers, data=postData)

		if not r.ok:
			r.raise_for_status()
			return None

		decoded = r.json()
		res = decoded # get whole list

	except Exception as e:
		print("Could not run VEP: ", e)
		return None
	d = []
	# return list of d
	for i in range(len(decoded)):
		d.append(get_output_list_pt2(res[i], hgvs[i], wt[i], mut[i], transcript[i]))
	return d

# temporary function contains rest of original for loop in get_output
def get_output_list_pt2(res, hgvs, wt="", mut="", transcript=""):
	# adjust for strand
	strand, indel_length, indel_bp = 1, 0, ""
	for key in res:
		# print(key, res[key])
		if (key == "strand" or "strand" in key) and (res[key] == -1):
			strand = -1
		if "del" in hgvs and (key == "allele_string" or "allele_string" in key):
			indel_info = res[key].strip().split("/")
			indel_length, indel_bp = len(indel_info[0]), str(indel_info[0])
		elif "ins" in hgvs and (key == "allele_string" or "allele_string" in key):
			indel_info = res[key].strip().split("/")
			indel_length, indel_bp = len(indel_info[1]), str(indel_info[1])
		elif "dup" in hgvs and (key == "allele_string" or "allele_string" in key):
			indel_info = res[key].strip().split("/")
			indel_length, indel_bp = len(indel_info[1]), str(indel_info[1])
	if strand == -1 and len(wt) != 0 and len(mut) != 0:
		splice_wt, splice_mut = nucleotide_change[wt], nucleotide_change[mut]
	else:
		splice_wt, splice_mut = wt, mut

	# find SIFT and Polyphen scores
	try:
		sift_score = [False]
		for item in res["transcript_consequences"]:
			if "sift_score" in item:
				sift_score = item["sift_score"]
				break
	except Exception as e:
		print("Could not retrieve sift_score: ", e)
		sift_score = [False]
	try:
		polyphen_score = [False]
		for item in res["transcript_consequences"]:
			if "polyphen_score" in item:
				polyphen_score = item["polyphen_score"]
				break
	except Exception as e:
		print("Could not retrieve polyphen_score: ", e)
		polyphen_score = [False]

	# find CADD score
	try:
		cadd_score = [False]
		for item in res["transcript_consequences"]:
			if "cadd_phred" in item:
				cadd_score = item["cadd_phred"]
				break
	except Exception as e:
		print("Could not retrieve cadd_score: ", e)
		cadd_score = [False]

	# find these scores (maxentscan, spliceAI)
	try:
		maxentscan_ref, ref_found = [False], False
		maxentscan_alt, alt_found = [False], False
		maxentscan_diff, diff_found = [False], False
		spliceai_pred, pred_found = [False], False
		for item in res["transcript_consequences"]:
			if "maxentscan_ref" in item and not ref_found:
				maxentscan_ref = item["maxentscan_ref"]
				ref_found = True
			if "maxentscan_alt" in item and not alt_found:
				maxentscan_alt = item["maxentscan_alt"]
				alt_found = True
			if "maxentscan_diff" in item and not diff_found:
				maxentscan_diff = item["maxentscan_diff"]
				diff_found = True
			if "spliceai_pred" in item and not pred_found:
				spliceai_pred = item["spliceai_pred"]
				pred_found = True
	except Exception as e:
		print("Could not retrieve some MES score: ", e)

	# look for start and chr of mutation
	try:
		gen_start = res["start"]
		gen_start = str(int(gen_start) - 1)  # Adjust for maxentscan perl script
		if "seq_region_name" in res:
			gen_chr = res["seq_region_name"]
		elif "colocated_variants" in res:
			gen_chr = res["colocated_variants"][0]["seq_region_name"]
		else:
			gen_chr = "un"
	except Exception as e:
		print("Could not retrieve some needed location information: ", e)

	# look for most severe consequence
	try:
		consequence = "Not found"
		if "most_severe_consequence" in res:
			consequence = res["most_severe_consequence"]
	except Exception as e:
		print("Could not retrieve most severe consequence: ", e)

	# look for amino acid change
	try:
		aa_change, codons, cdna = "Not found", "not found", "Not found"
		for item in res["transcript_consequences"]:
			if "amino_acids" in item:
				aa_change = item["amino_acids"]
			if "codons" in item:
				codons = item["codons"]
			if "cdna_start" in item:
				cdna = item["cdna_start"]
	except Exception as e:
		print("Could not retrieve amino acid change: ", e)

	# create gnomad variant string
	gnomad_variant = ""
	if gen_chr and gen_start and len(splice_wt) != 0 and len(splice_mut) != 0:
		gnomad_variant = "{}-{}-{}-{}".format(gen_chr, str(int(gen_start) + 1), splice_wt, splice_mut)

	# # get spliceai prediction from spliceai-lookup and mes scores
	spliceai_pred = spliceai(gnomad_variant)
	mes_dict_5 = mes5_runner(gen_start, gen_chr, mut, strand, hgvs, indel_length, indel_bp)
	mes_dict_3 = mes3_runner(gen_start, gen_chr, mut, strand, hgvs, indel_length, indel_bp)

	# Determine if this mutation is amenable
	coding_impact, splicing_impact, consequence, mes_analyses, spliceai_analyses = determine_impacts(gen_chr, gen_start,
																									 hgvs, consequence,
																									 sift_score,
																									 polyphen_score,
																									 maxentscan_ref,
																									 maxentscan_alt,
																									 maxentscan_diff,
																									 spliceai_pred,
																									 mes_dict_5,
																									 mes_dict_3, strand)

	# IMPT:
	# mes_analyses is a dictionary that maps key: "key,ref,alt,delta" to a tuple value: (x, y, z) where x=high/moderate/low prob, y=donor/acceptor, z=weakened/strengthened
	# spliceai_analyses is a dictionary that maps key: "pos_delta,prob_delta" to a tuple value: (x, y), where x=acceptor/donor, y=loss/gain
	print(mes_analyses, spliceai_analyses)
	aso_prediction = determine_amenability(mes_analyses, spliceai_analyses, hgvs, consequence, sift_score,
										   polyphen_score, strand)

	d = {
		"aso_prediction": aso_prediction,
		"coding_impact": coding_impact,
		"splicing_impact": splicing_impact,
		"sift_score": sift_score,
		"polyphen_score": polyphen_score,
		"maxentscan_ref": maxentscan_ref,
		"maxentscan_alt": maxentscan_alt,
		"maxentscan_diff": maxentscan_diff,
		"spliceai_pred": spliceai_pred,
		"gnomad_variant": gnomad_variant,
		"cadd_score": cadd_score,
		"strand": str(strand),
		"aa_change": aa_change,
		"codons": codons,
		"cdna": cdna
	}
	return d


def main():
	# tup2 = get_output("NPC1:c.1554-1009G>A", "G", "A")
	# tup2 = get_output("LMNA:c.640-10A>G", "A", "G")
	# tup2 = get_output("LMNA:c.1824C>T", "C", "T")
	# tup2 = get_output("APC:c.532-1000delGT", "", "")
	# tup2 = get_output("LMNA:c.1488+5G>C", "G", "C")
	# tup2 = get_output("DGAT1:c.1013_1015delTCT")
	# tup2 = get_output("MUSK:c.1189_1190insCGT")
	# tup2 = get_output("WDR81:c.2292_2309del18")
	# tup2 = get_output("MGP:c.94+1G>A", "G", "A")
	# tup2 = get_output("ABCA4:c.5461-10T>C", "T", "C")
	# tup2 = get_output("B4GALNT1:c.532-1G>C", "G", "C")
	# tup2 = get_output("ALS2:c.2581-1G>C", "G", "C")
	tup2 = get_output("ATM:c.8500T>A", "T", "A")

	print(tup2)

	# variants = [
	# ("ATM:c.7865C>T", "C", "T"),
	# ("ATM:c.6203T>C", "T", "C"),
	# ("ATM:c.5763-1050A>G", "A", "G"),
	# ("ATM:c.3994-159A>G", "A", "G"),
	# ("ATM:c.3489C>T", "C", "T"),
	# ("ATM:c.4801A>G", "A", "G"),
	# ("ATM:c.7318A>G", "A", "G"),
	# ("ATM:c.6348-987G>C", "G", "C"),
	# ("ATM:c.7914G>T", "G", "T"),
	# ("ATM:c.2639-21A>G", "A", "G"),
	# ("ATM:c.3529T>C", "T", "C"),
	# ("ATM:c.6573-15T>G", "T", "G"),
	# ("ATM:c.8494C>T", "C", "T"),
	# ("ATM:c.5228C>T", "C", "T"),
	# ("ATM:c.6047A>G", "A", "G"),
	# ("ATM:c.8495G>T", "G", "T"),
	# ("ATM:c.8486C>T", "C", "T"),
	# ("ATM:c.743G>T", "G", "T"),
	# ("ATM:c.496+5G>A", "G", "A"),
	# ("ATM:c.2250G>A", "G", "A")
	# ]

	# output_file = open("atm_variants.tsv", "w+")
	# output_file.write("variant\taso_prediction\tcoding_impact\tsplicing_impact\n")

	# for variant in variants:
	# 	print(variant)
	# 	tup = get_output(variant[0], variant[1], variant[2])
	# 	output_file.write("{}\t{}\t{}\t{}\n".format(variant[0], tup["aso_prediction"], tup["coding_impact"], tup["splicing_impact"]))
	# output_file.close()

	# d = {
	# 	"aso_prediction": aso_prediction,
	# 	"coding_impact": coding_impact,
	# 	"splicing_impact": splicing_impact,
	# 	"sift_score": sift_score,
	# 	"polyphen_score": polyphen_score,
	# 	"maxentscan_ref": maxentscan_ref,
	# 	"maxentscan_alt": maxentscan_alt,
	# 	"maxentscan_diff": maxentscan_diff,
	# 	"spliceai_pred": spliceai_pred,
	# 	"gnomad_variant": gnomad_variant,
	# 	"cadd_score": cadd_score,
	# 	"strand": str(strand),
	# 	"aa_change": aa_change,
	# 	"codons": codons,
	# 	"cdna": cdna
	# }

if __name__== "__main__":
	main()

