import os
import subprocess
from multiprocessing import Pool

# ================= Configuration =================
MAX_PARALLEL_JOBS = 6         # 同时跑的任务数 (方案 B)
CORES_PER_JOB = 1              # 每个任务给 2 核
ENERGY_FILE = "energy_all.txt" # 你实际输出能量的文件名
EXECUTABLE = "lmp_mpinep_new"         # LAMMPS 执行文件名
# =================================================

def get_folders():
    folders = [d for d in os.listdir() if d.startswith("iter00_") and os.path.isdir(d)]
    folders.sort(key=lambda x: int(x.split('_')[1]))
    return folders

def run_lammps_task(dir_name):
    """
    单个文件夹的执行逻辑
    """
    dir_path = os.path.join(os.getcwd(), dir_name)
    output_path = os.path.join(dir_path, ENERGY_FILE)

    # 1. 检查是否已算过
    if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
        return f"[Skip] {dir_name} already finished."

    # 2. 构造运行命令 (mpirun -n 2 lmp_mpi < in.GSFE1)
    # 使用 -in 更加稳健，并将输出重定向
    cmd = f"mpirun -n {CORES_PER_JOB} {EXECUTABLE} -in in.elastic"
    
    try:
        # 在子目录中执行
        with open(os.path.join(dir_path, "lammps.out"), "w") as out_f:
            subprocess.run(cmd, shell=True, cwd=dir_path, stdout=out_f, stderr=subprocess.STDOUT)
        return f"[Done] {dir_name} finished."
    except Exception as e:
        return f"[Error] {dir_name}: {str(e)}"

def main():
    all_dirs = get_folders()
    total = len(all_dirs)
    print(f"Total folders: {total}. Max parallel jobs: {MAX_PARALLEL_JOBS}")
    print("Starting local parallel execution...")

    # 创建进程池
    with Pool(processes=MAX_PARALLEL_JOBS) as pool:
        # imap 会返回一个迭代器，我们可以实时打印进度
        for result in pool.imap(run_lammps_task, all_dirs):
            print(result)

if __name__ == "__main__":
    main()
