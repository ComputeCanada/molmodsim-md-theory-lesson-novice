import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots()
fig.subplots_adjust(left=0.26)
fig.set_size_inches(9,11)
# Format: software/resources, parallel speed (ns/day), color, [serial speed (ns/day)]
benchmark=[
        [ 'AMBER/18, pmemd.MPI\n 32 tasks x 1 core',                  2.16, 'tab:green' ],        
	[ 'NAMD-multicore/2.14\n 1 task x 32 cores',          2.18, 'tab:red' ],
        [ 'AMBER/18, pmemd.MPI\n 64 tasks x 1 core',                  3.62, 'tab:green' ],
	[ 'NAMD-UCX/2.14\n 80 tasks x 1 core',                        4.67, 'tab:red' ],
        [ 'AMBER/18, pmemd.MPI\n128 tasks x 1 core',                  5.77, 'tab:green' ], 
	[ 'NAMD-multicore-cuda/2.14\n1 task x 8 cores, 1 GPU',        5.92, 'tab:red' ],
        [ 'AMBER/18, pmemd.MPI\n256 tasks x 1 core',                  8.78, 'tab:green' ], 
        [ 'NAMD-multicore-cuda/2.14\n1 task x 16 cores, 2 GPUs',      8.85, 'tab:red' ],
        [ 'NAMD-UCX/2.14\n160 tasks x 1 core',                        9.06, 'tab:red' ],
        [ 'GROMACS/2020.4, AVX2\n32 tasks x 1 core',                 10.27, 'tab:blue' ],
        [ 'GROMACS/2020.4, AVX512\n32 tasks x 1 core',               10.99, 'tab:blue' ],
        [ 'NAMD-multicore-cuda/2.14\n1 task x 40 cores, 2 GPUs',     11.17, 'tab:red' ],
        [ 'GROMACS/2020.4, CUDA\n2 tasks x 8 cores, 2 GPUs', 12.39, 'tab:blue' ],    
        [ 'GROMACS/2020.4, AVX2\n64 tasks x 1 core',                 17.03, 'tab:blue' ],
        [ 'GROMACS/2020.4, AVX512\n64 tasks x 1 core',               20.76, 'tab:blue' ],
        [ 'GROMACS/2020.4, CUDA\n1 task x 8 cores, 1 GPU',           22.25, 'tab:blue' ],       
        [ 'NAMD/3.0_alpha8, CUDA\n1 task x 1 cores, 1 GPU',          23.10, 'tab:red' ], 
        [ 'GROMACS/2020.4, AVX512\n32 task x 4 cores',               27.27, 'tab:blue' ],
        [ 'GROMACS/2020.4, AVX2\n128 tasks x 1 core',                33.32, 'tab:blue' ],
        [ 'GROMACS/2020.4, AVX512\n128 tasks x 1 core',              36.56, 'tab:blue' ],
        [ 'AMBER/18, pmemd.cuda.MPI\n2 tasks x 1 core, 2 GPUs',      36.88, 'tab:green' ],
        [ 'AMBER/18, pmemd.cuda\n1 task x 1 core, 1 GPU',            41.77, 'tab:green' ],
        [ 'GROMACS/2020.4, AVX2\n256 tasks x 1 core',                52.10, 'tab:blue' ],        
        [ 'GROMACS/2020.4, AVX512\n256 tasks x 1 core',              56.10, 'tab:blue' ],
]

software=[s[0] for s in benchmark]
rate=[s[1] for s in benchmark]
clr=[s[2] for s in benchmark]
gmx_serial=0.522


y_pos = np.arange(len(rate))
ax.barh(y_pos, rate, color=clr)
ax.set_yticks(y_pos)
ax.set_yticklabels(software)
# ax.invert_yaxis() 
ax.set_xlabel('Performance, ns/day')
ax.set_title('NPT, 1 fs, 300,000 atoms\nCPU: Xeon Gold 6248 @ 2.50GHz, GPU: Tesla V100')
ax.minorticks_on()
ax.grid(which='major', linestyle='-', linewidth='0.7', alpha=0.6, axis='x')
ax.grid(which='minor', linestyle='-', linewidth='0.5', alpha=0.2, axis='x')
ax.set_axisbelow(True)


i=10
ax.text(rate[i]-6, i-0.1, "{:.1f}".format(rate[i]*100/(gmx_serial*32))+"%", color="white")
i=14
ax.text(rate[i]-6, i-0.1, "{:.1f}".format(rate[i]*100/(gmx_serial*64))+"%", color="white")
i=19
ax.text(rate[i]-6, i-0.1, "{:.1f}".format(rate[i]*100/(gmx_serial*128))+"%", color="white")
i=23
ax.text(rate[i]-6, i-0.1, "{:.1f}".format(rate[i]*100/(gmx_serial*256))+"%", color="white")


plt.savefig('MD-benchmarks.svg', dpi=600)

#plt.show()

