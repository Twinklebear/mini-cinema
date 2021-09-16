#!/bin/bash

get_physical_cores() {
    echo `grep "^cpu\\scores" /proc/cpuinfo | uniq | awk '{print $4}'`
}

get_logical_cores() {
    echo `grep "^processor" /proc/cpuinfo | uniq | wc -l`
}

set_ospray_env_vars() {
    #export I_MPI_PIN_RESPECT_CPUSET=0
    #export I_MPI_PIN_RESPECT_HCA=0
    export I_MPI_PIN_DOMAIN=omp
    #export OSPRAY_SET_AFFINITY=0

    CPU_MODEL=`cat /proc/cpuinfo | grep "model name" | uniq`
    if [ -z "$OSPRAY_NUM_THREADS" ]; then
        # If we're on Xeon use all logical cores, if on Phi use only phyiscal cores
        if [ -n "`echo $CPU_MODEL | grep 'Xeon(R)'`" ]; then
            echo "CPU $CPU_MODEL is Xeon, using all logical cores"
            export OSPRAY_NUM_THREADS=$(get_logical_cores)
            export I_MPI_PIN_PROCESSOR_LIST=all
        else
            echo "CPU $CPU_MODEL is Xeon Phi, using all physical cores"
            export OSPRAY_NUM_THREADS=$(get_physical_cores)
            export I_MPI_PIN_PROCESSOR_LIST=allcores
        fi
        echo "I_MPI_PIN_PROCESSOR_LIST = $I_MPI_PIN_PROCESSOR_LIST"
        echo "OSPRAY_NUM_THREADS = $OSPRAY_NUM_THREADS"
    else
        echo "OSPRAY_NUM_THREADS set by user to $OSPRAY_NUM_THREADS"
    fi

    export OMP_NUM_THREADS=$OSPRAY_NUM_THREADS
    export MPICH_MAX_THREAD_SAFETY=multiple
}

set_ospray_env_vars

