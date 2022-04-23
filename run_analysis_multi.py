import os
import concurrent.futures
import subprocess

def main():
    run_info = [line.rstrip() for line in open('run_analysis.txt')]

    num_iterations = []
    run_strings = []

    total_threads = 0

    for i in range(len(run_info)):
        if i % 2:
            run_strings.append(run_info[i])
        else:
            num_iterations.append(int(run_info[i]))

    to_run = []


    for i in range(len(num_iterations)):
        # find number of digits to be presented to the
        total_digits = num_iterations[i] // 10

        total_threads += num_iterations[i]

        to_run.extend([run_strings[i]+f'-{n}'.zfill(total_digits) +f' --seed={n}' for n in range(num_iterations[i])])

    for j in to_run:
        print(j)
    # print(to)
    # exit()
    print(f'{total_threads} threads to spawn')

    # print(to_run)
    # exit()

    i = 0
    # i = 0
    # for f in to_run:
    #    p = subprocess.Popen(["/c/Users/DNA_is_cool/Documents/oxDNA/build/bin/oxDNA",f])
    #    i=1+i
    #    print(f'file number {i+1} has been run')

    with concurrent.futures.ThreadPoolExecutor(max_workers=total_threads) as executor:
        results = [executor.submit(run_analysis_script, runstring) for runstring in to_run]
            # results = executor.starmap(rollout_omt,[(G,d_e,node) for d_e in to_remove])
        for f in concurrent.futures.as_completed(results):
            try:
                c = f.result()
                i += 1
                print(f'{i} simulations are complete')
            except (KeyError):
                print(f'Error')
            except Exception as e:
                print(e.output)
            
    print(f'All simulations are complete!')

def run_analysis_script(runstring):
    p = subprocess.run(runstring, shell=True, stdout=subprocess.DEVNULL).returncode
    return p

if __name__ == '__main__':
    main()


