import threading

t1 = threading.Thread(target=get_output, args=(hgvs, wt, mut,))

t1.start()
