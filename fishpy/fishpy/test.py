import time
time_start = time.time()
n = 10
i = 1
for x in range(n):
    time.sleep(2.4)
    time_tmp = time.time()
    ETA = (n - i)*(time_tmp - time_start) / i
    print("{:3.1f}% done, ETA : {:2.0f} min {:2.0f} secs".format(i/n * 100, ETA // 60, ETA % 60), end = "\r" )
    i +=1
