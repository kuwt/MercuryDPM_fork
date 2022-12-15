import multiprocessing
import psutil

numFreeCores = int(multiprocessing.cpu_count()*(0.5-0.01*psutil.cpu_percent(interval=1)))
print(numFreeCores)
