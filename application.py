import os
import asyncio

from flask import Flask, flash, jsonify, redirect, render_template, request, session, send_file, Response, make_response
import os
import pandas as pd
from flask_session import Session
from tempfile import mkdtemp
from werkzeug.exceptions import default_exceptions, HTTPException, InternalServerError
from werkzeug.security import check_password_hash, generate_password_hash
from spliceCheck_new2 import *
from helpers import apology
from datetime import datetime
from gevent import monkey


# Configure application
application = app = Flask(__name__)
app.secret_key = '84e1a45ba1acb48d7d7c988329bcf07ab303b3a85e68889cc54b7cbd049b272922d3857f8deebfcd1ea7'

# Ensure templates are auto-reloaded
app.config["TEMPLATES_AUTO_RELOAD"] = True


# Ensure responses aren't cached

@app.after_request
def after_request(response):
    response.headers["Cache-Control"] = "no-cache, no-store, must-revalidate"
    response.headers["Expires"] = 0
    response.headers["Pragma"] = "no-cache"
    return response


# Configure session to use filesystem (instead of signed cookies)
# app.config["SESSION_FILE_DIR"] = mkdtemp()
# app.config["SESSION_PERMANENT"] = False
# app.config["SESSION_TYPE"] = "filesystem"
# Session(app)

@app.route("/")
def index():
    return render_template("index.html")


@app.route("/list")
def list_page():
    return render_template("listpage.html")


@app.route("/about")
def about():
    """Show information about the tool"""
    # Go to the page
    return render_template("about.html")


@app.route("/output", methods=["POST"])
def getOutput():
    """Call tools to get output"""

    # Get user input from form
    gene = request.form.get("gene")
    cdna = request.form.get("position")
    wt = request.form.get("wt")
    mut = request.form.get("mut")
    transcript = request.form.get("transcript")

    # Check if it's a possible gene name
    if not gene.isalnum():
        return apology("Try fixing the gene name")

    # Check mutation

    if "del" in cdna or "ins" in cdna:
        wt, mut = "", ""
    elif str(wt) == str(mut):
        return apology("Double check that the mutation is correct")

    # Check if the coordinate format is correct
    cdna_check = list(cdna)
    nums = ''.join(cdna_check[2:])
    if cdna_check[0] != 'c' or cdna_check[1] != '.' or len(nums) < 1:
        return apology("Make sure the coordinate format is correct")

    # HGVS notation that can be used to call VEP (gene:c.123wt>mut)
    # hgvs = gene + ':' + cdna + wt + '>' + mut
    if len(wt) == 0 and len(mut) == 0:
        hgvs = "{}:{}".format(gene.upper(), cdna)
    else:
        hgvs = "{}:{}{}>{}".format(gene.upper(), cdna, wt, mut)

    # # THIS IS WHERE GET_OUTPUT IN spliceCheck_new.py CAN BE CALLED#####################
    res = get_output(hgvs, wt, mut, transcript)
    if not res:
        return apology("I could not check to see if the mutation is amenable. Please wait ")
    elif res == "1":
        return apology("Could not run VEP. Please wait ")
    elif res == "2":
        return apology("Could not run VEP. Please wait ")
    for key in res:  # change [False] to "None found"
        if type(res[key]) is list:
            res[key] = "None found."
    if res["strand"] == "1":
        res["strand"] = "Positive"
    else:
        res["strand"] = "Negative"
    scores = [
        "Variant: " + hgvs,
        "Amino acid change: {} ({})".format(str(res["aa_change"]), str(res["codons"])),
        "cDNA coordinate: " + str(res["cdna"]),
        "SIFT: " + str(res["sift_score"]),
        "PolyPhen: " + str(res["polyphen_score"]),
        "CADD Phred: " + str(res["cadd_score"]),
        "VEP MaxEntScan Ref: " + str(res["maxentscan_ref"]),
        "VEP MaxEntScan Alt: " + str(res["maxentscan_alt"]),
        "VEP MaxEntScan Diff: " + str(res["maxentscan_diff"]),
        "SpliceAI: " + str(res["spliceai_pred"]),
        "Strand: " + str(res["strand"]),
        "Transcript ID: " + str(transcript)
    ]
    if "stop-gain" in res["coding_impact"].lower():
        ans = res["splicing_impact"].strip()
    else:
        ans = "{} {}".format(res["coding_impact"], res["splicing_impact"]).strip()
    ans = ans.split(". ")  # should separate list into parts to bullet - make sure this is robust
    print(ans)

    # Call literature search
    literature = get_literature(hgvs)
    if not literature:  # double check returnable result
        literature = []

    # Call clinvar search
    clinvar_id = get_clinvar(hgvs)
    if len(clinvar_id) == 0:  # double check returnable result
        clinvar_id = ""

    return render_template("output.html", hgvs=hgvs, literature=literature, clinvar=clinvar_id,
                           pred=res["aso_prediction"], ans=ans, vep_output=scores, gene_id=gene.upper(),
                           var=res["gnomad_variant"], wt5="None", wt3="None", mut5="None", mut3="None")


# wt5=res["wt5mes"], wt3=res["wt3mes"], mut5=res["mut5mes"], mut3=res["mut3mes"])

# def getOutputList():
#     # Get user input from form
#     file = request.files["file"]
#     file.seek(0)
#     contents = file.read().decode("utf-8")
#     variant_list = contents.strip().split("\n")
#     scores = []
#     for hgvs in variant_list:
#         res = "0"
#
#         # input checking
#         hgvs_comps = hgvs.strip().split(">")
#         if len(hgvs_comps) == 2:
#             wt, mut = hgvs_comps[0][-1], hgvs_comps[1].strip()
#             cdna = hgvs_comps[0][0:-1]
#             if len(cdna.strip().split(":")) != 2:
#                 res = "1"
#             elif cdna.strip().split(".")[0][-1] != "c":
#                 res = "1"
#             if str(wt) == str(mut):
#                 res = "1"
#         elif len(hgvs_comps) == 1:   # insertions/deletions/duplications - frameshift
#             wt, mut = "", ""
#             if len(hgvs_comps[0].strip().split(":")) != 2:
#                 res = "1"
#             elif hgvs_comps[0].strip().split(".")[0][-1] != "c":
#                 res = "1"
#             elif "dup" not in hgvs and "del" not in hgvs and "ins" not in hgvs:
#                 res = "1"
#         else:
#             res = "1"
#
#         if res == "1":
#             scores.append((hgvs, "Please double check that the mutation was inputted correctly."))
#         else:
#             res = get_output(hgvs, wt, mut)
#             if not res:
#                 scores.append((hgvs, "Could not run VEP on this variant. Please double check that the mutation was inputted correctly."))
#             else:
#                 for key in res: # change [False] to "None found"
#                     if type(res[key]) is list:
#                         res[key] = "None found."
#                 str_result = "coding_impact: {} | splicing_impact: {} | SIFT: {} | PolyPhen: {} | VEP MaxEntScan Ref: {} | VEP MaxEntScan Alt: {} | VEP MaxEntScan Diff: {} | SpliceAI: {} | gnomad_variant: {} ".format(res["coding_impact"], res["splicing_impact"], res["sift_score"], res["polyphen_score"], res["maxentscan_ref"], res["maxentscan_alt"], res["maxentscan_diff"], res["spliceai_pred"], res["gnomad_variant"])
#                 scores.append((hgvs, str_result))
#     return render_template("outputFile.html", vep_output=scores)

# need to fix
@app.route("/outputFile", methods=["POST"])
def getOutputList():
    # Get user input from form
    file = request.files["file"]
    file.seek(0)
    contents = file.read().decode("utf-8")
    # separate variant by lines and use transcript
    contents = contents.splitlines()
    variant_list = []
    transcript_list = []
    for content in contents:
        variant = content.split(" ")[0]
        if len(variant) == len(content):
            transcript = ""
        else:
            transcript = content.split(" ")[1]
        variant_list.append(variant)
        transcript_list.append(transcript)

    # initialize list of results to value "0"
    res = ["0" for _ in range(len(variant_list))]
    # empty scores list
    scores = []
    # empty input lists (match to input of get_output
    hgvs_list = []
    wt_list = []
    mut_list = []

    # input checking for loop
    for i, hgvs in enumerate(variant_list):

        # input checking
        hgvs_comps = hgvs.strip().split(">")  # check for substitution
        if len(hgvs_comps) == 2:
            wt, mut = hgvs_comps[0][-1], hgvs_comps[1].strip()
            cdna = hgvs_comps[0][0:-1]
            if len(cdna.strip().split(":")) != 2:
                res[i] = "1"
            elif cdna.strip().split(".")[0][-1] != "c":
                res[i] = "1"
            if str(wt) == str(mut):
                res[i] = "1"
        elif len(hgvs_comps) == 1:  # insertions/deletions/duplications - frameshift
            wt, mut = "", ""
            if len(hgvs_comps[0].strip().split(":")) != 2:
                res[i] = "1"
            elif hgvs_comps[0].strip().split(".")[0][-1] != "c":
                res[i] = "1"
            elif "dup" not in hgvs and "del" not in hgvs and "ins" not in hgvs:
                res[i] = "1"
        else:
            res[i] = "1"
        # if it fails add failure text to scores output: # TODO: it will out of order, if it matters
        if res[i] == "1":
            scores.append((hgvs, "Please double check that the mutation was inputted correctly."))
        # else add to splicecheck_input
        else:
            hgvs_list.append(hgvs)
            wt_list.append(wt)
            mut_list.append(mut)

    # for all input ready to execute, run get_output_list
    output_res = get_output_list(hgvs_list, wt_list, mut_list, transcript_list)
    # todo: check that output_res not none
    for i, result in enumerate(output_res):
        if not result:
            scores.append((hgvs_list[i],
                           "Could not run VEP on this variant. Please double check that the mutation was inputted correctly."))
        else:
            for key in result:  # change [False] to "None found"
                if type(result[key]) is list:
                    result[key] = "None found."
            str_result = "coding_impact: {} | splicing_impact: {} | SIFT: {} | PolyPhen: {} | VEP MaxEntScan Ref: {} | VEP MaxEntScan Alt: {} | VEP MaxEntScan Diff: {} | SpliceAI: {} | gnomad_variant: {} ".format(
                result["coding_impact"], result["splicing_impact"], result["sift_score"], result["polyphen_score"],
                result["maxentscan_ref"], result["maxentscan_alt"], result["maxentscan_diff"], result["spliceai_pred"],
                result["gnomad_variant"])
            scores.append((hgvs_list[i], str_result))
    session['scores'] = scores
    # todo: log
    # df = format_csv(scores)

    return render_template("outputFile.html", vep_output=scores)

@app.route('/export.csv')
def download_csv():
    scores = session['scores']
    df = format_csv(scores)
    print(df)
    return Response(
        # df.to_csv(index=False, header=False),
        df.to_csv(index=False, header=True),
        mimetype="text/csv",
        headers={"Content-disposition":
                     "attachment; filename=export.csv"})

# download logs of last queries
@app.route('/log.csv')
def download_log():
    # just return the log
    return Response("output/log.csv", mimetype="text/csv", headers={"Content-disposition": "attachment; filename=export.csv"})

# helper function to make csv from list of lists
def format_csv(scores):
    # Returns a datetime object containing the local date and time
    dateTimeObj = datetime.now()
    timestampStr = dateTimeObj.strftime("%d-%b-%Y (%H:%M:%S)")
    # reformat list of lists
    new_scores = []
    for tup in scores:
        # remove the carriage returns...
        # row = [tup[0].strip()] + [x.strip() for x in tup[1].split('|')] + [timestampStr]
        row = [tup[0].strip()] + [x.split(":")[1] for x in tup[1].split(' | ')] + [timestampStr]
        new_scores.append(row)
    df = pd.DataFrame(new_scores, columns=['Variant', 'Coding Impact', 'Splicing Impact', 'SIFT', 'Polyphen', 'VEP MaxEntScan Ref', 'VEP MaxEntScan Alt','VEP MaxEntScan Diff','SpliceAI','gnomad_variant','Time'])
    return df

if __name__ == "__main__":
    # context = ('server.crt', 'privatekey.pem')  # certificate and key files
    monkey.patch_all()
    app.run(debug=True)

    # app.run(debug=True, port=443)
    # app.run(debug=True, ssl_context=context)
    # https://stackoverflow.com/questions/47108504/flask-default-port-number
    # app.run(host='0.0.0.0', port=80)

# # Run VEP on the variant
# v = os.popen(' ./tools/ensembl-vep/vep --database -id "%s" -o STDOUT --force --sift s --polyphen s' % hgvs, 'r')
# # Read in the output
# vep = v.readlines()
# # Make sure this worked
# if len(vep) < 36:
#     return apology("VEP failed to run, please re-enter your variant")
# # The output is long - we only care about the 35th line
# vep = ''.join(vep[35])
# # Split this line up into a list
# vep = vep.split('\t')

# # The second item in the list is the genomic coordinates chr:number
# gcoords = vep[1]
# # The last item is a list of scores we care about
# scores = vep[-1].split(';')
# # The last item in this list is Polyphen -- we just want the score
# polyP = scores[-1].strip()
# polyP = polyP.split('=')
# polyP = float(polyP[-1])
# # The second to last is SIFT -- we just want the score
# siftScore = scores[-2]
# siftScore = siftScore.split('=')
# siftScore = siftScore[-1]
# # gnomAD requires the variant in a specific genomic coordinate form
# gncoords = gcoords.split(':')
# variant = gncoords[0]+'-'+gncoords[1]+'-'+wt+'-'+mut

# # MAX ENT SCAN
# # Run 5' Max Ent Scan
# # Get the 9 bp sequence around the mutation
# me5low = int(gncoords[1])-4
# me5high = int(gncoords[1])+5
# me5range = 'chr'+gncoords[0]+":"+str(me5low)+'-'+str(me5high)

# # This just pulls the sequence from an online database
# f = os.popen('./tools/twoBitToFa http://hgdownload.cse.ucsc.edu/gbdb/hg19/hg19.2bit:%s stdout' % me5range, 'r')
# # Get the output
# lst = f.readlines()

# # Make sure this worked
# if len(lst) < 2:
#     return apology("Something went wrong. Please re-enter your coordinates %s" % (''.join(lst)))

# # The second value is the sequence
# seq9 = (lst[1].strip()).lower()

# # Call MaxEntScan 5' on this sequence
# wt5mes = mes5(seq9)

# # Change the middle value to the mutation
# seq9m = seq9[0:4]+mut+seq9[6:9]
# # Call MaxEntScan 5' on this mutant seq
# mut5mes = mes5(seq9)

# # Run 3' Max Ent Scan
# # Get the  bp sequence around the mutation
# me23low = int(gncoords[1])-11
# me23high = int(gncoords[1])+12
# me23range = 'chr'+gncoords[0]+":"+str(me23low)+'-'+str(me23high)

# # This just pulls the sequence from an online database
# f = os.popen('./tools/twoBitToFa http://hgdownload.cse.ucsc.edu/gbdb/hg19/hg19.2bit:%s stdout' % me23range, 'r')
# # Get the output
# lst2 = f.readlines()

# # Make sure this worked
# if len(lst2) < 2:
#     return apology("Something went wrong. Please re-enter your coordinates%s" % (''.join(lst2)))

# # The second value is the sequence
# seq23 = lst2[1].strip().lower()

# # Call MaxEntScan 5' on this sequence
# wt3mes = mes3(seq23)

# # Change the middle value to the mutation
# seq23m = seq23[0:11]+mut+seq23[13:23]
# # Call MaxEntScan 5' on this mutant seq
# mut3mes = mes3(seq23)

# # Determine if this mutation is ammenable
# answer = fixable(siftScore, polyP, wt5mes, mut5mes, wt3mes, mut3mes)
