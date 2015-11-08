main:
	nvcc -arch=compute_30 -code=sm_30 ./main.cu -o main

run: main
	srun -J hcd -o hcd.o%j -p gpudev -n 1 -t 00:10:00 ./main

clean:
	rm -rf main

.PHONY: run clean
