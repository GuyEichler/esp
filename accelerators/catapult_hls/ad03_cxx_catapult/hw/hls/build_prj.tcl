# Copyright (c) 2011-2020 Columbia University, System Level Design Group
# SPDX-License-Identifier: Apache-2.0

#
# Accelerator
#

set ACCELERATOR "ad03_cxx_catapult"
set PLM_HEIGHT 128
set PLM_WIDTH 32
set PLM_SIZE [expr ${PLM_WIDTH}*${PLM_HEIGHT}]

set uarch "basic"

#
# Technology-dependend reports and project dirs.
#

file delete -force -- $ACCELERATOR\_$uarch\_fx$PLM_WIDTH\_dma$DMA_WIDTH
file delete -force -- $ACCELERATOR\_$uarch\_fx$PLM_WIDTH\_dma$DMA_WIDTH.css
project new -name $ACCELERATOR\_$uarch\_fx$PLM_WIDTH\_dma$DMA_WIDTH
set CSIM_RESULTS "./tb_data/catapult_csim_results.log"
set RTL_COSIM_RESULTS "./tb_data/catapult_rtl_cosim_results.log"

#
# Reset the options to the factory defaults
#

solution new -state initial
solution options defaults

solution options set Flows/ModelSim/VLOG_OPTS {-suppress 12110}
solution options set Flows/ModelSim/VSIM_OPTS {-t ps -suppress 12110}
solution options set Flows/DesignCompiler/OutNetlistFormat verilog
solution options set /Input/CppStandard c++11
#solution options set /Input/TargetPlatform x86_64

set CATAPULT_VERSION  [string map { / - } [string map { . - } [application get /SYSTEM/RELEASE_VERSION]]]
solution options set Cache/UserCacheHome "catapult_cache_$CATAPULT_VERSION"
solution options set Cache/DefaultCacheHomeEnabled false
solution options set /Flows/SCVerify/DISABLE_EMPTY_INPUTS true

flow package require /SCVerify

flow package option set /SCVerify/USE_CCS_BLOCK true
flow package option set /SCVerify/USE_QUESTASIM true
flow package option set /SCVerify/USE_NCSIM false

#options set Flows/OSCI/GCOV true
#flow package require /CCOV
#flow package require /SLEC
#flow package require /CDesignChecker

#directive set -DESIGN_GOAL area
##directive set -OLD_SCHED false
#directive set -SPECULATE true
#directive set -MERGEABLE true
directive set -REGISTER_THRESHOLD 8192
#directive set -MEM_MAP_THRESHOLD 32
#directive set -LOGIC_OPT false
#directive set -FSM_ENCODING none
#directive set -FSM_BINARY_ENCODING_THRESHOLD 64
#directive set -REG_MAX_FANOUT 0
#directive set -NO_X_ASSIGNMENTS true
#directive set -SAFE_FSM false
#directive set -REGISTER_SHARING_MAX_WIDTH_DIFFERENCE 8
#directive set -REGISTER_SHARING_LIMIT 0
#directive set -ASSIGN_OVERHEAD 0
#directive set -TIMING_CHECKS true
#directive set -MUXPATH true
#directive set -REALLOC true
#directive set -UNROLL no
#directive set -IO_MODE super
#directive set -CHAN_IO_PROTOCOL standard
#directive set -ARRAY_SIZE 1024
#directive set -REGISTER_IDLE_SIGNAL false
#directive set -IDLE_SIGNAL {}
#directive set -STALL_FLAG false
##############################directive set -TRANSACTION_DONE_SIGNAL true
#directive set -DONE_FLAG {}
#directive set -READY_FLAG {}
#directive set -START_FLAG {}
#directive set -BLOCK_SYNC none
#directive set -TRANSACTION_SYNC ready
#directive set -DATA_SYNC none
#directive set -CLOCKS {clk {-CLOCK_PERIOD 0.0 -CLOCK_EDGE rising -CLOCK_UNCERTAINTY 0.0 -RESET_SYNC_NAME rst -RESET_ASYNC_NAME arst_n -RESET_KIND sync -RESET_SYNC_ACTIVE high -RESET_ASYNC_ACTIVE low -ENABLE_ACTIVE high}}
directive set -RESET_CLEARS_ALL_REGS true
#directive set -CLOCK_OVERHEAD 20.000000
#directive set -OPT_CONST_MULTS use_library
#directive set -CHARACTERIZE_ROM false
#directive set -PROTOTYPE_ROM true
#directive set -ROM_THRESHOLD 64
#directive set -CLUSTER_ADDTREE_IN_COUNT_THRESHOLD 0
#directive set -CLUSTER_OPT_CONSTANT_INPUTS true
#directive set -CLUSTER_RTL_SYN false
#directive set -CLUSTER_FAST_MODE false
#directive set -CLUSTER_TYPE combinational
#directive set -COMPGRADE fast

#set CLOCK_PERIOD 12.5

# Flag to indicate SCVerify readiness
set can_simulate 1

# Design specific options.

#
# Flags
#

solution options set Flows/QuestaSIM/SCCOM_OPTS {-64 -g -x c++ -Wall -Wno-unused-label -Wno-unknown-pragmas}

#
# Input
#

solution options set /Input/SearchPath { \
    ../tb \
    ../inc \
    ../firmware \
    ../firmware/ap_types \
    ../firmware/nnet_utils \
    ../firmware/weights \
    ../src \
    ../../../common/inc }

# Add source files.
solution file add ../src/ad03_cxx_catapult.cpp -type C++
solution file add ../inc/ad03_cxx_catapult.hpp -type C++
solution file add ../firmware/anomaly_detector.cpp -type C++
solution file add ../firmware/anomaly_detector.h -type C++
solution file add ../tb/main.cpp -type C++ -exclude true

solution file set ../src/ad03_cxx_catapult.cpp -args -DDMA_WIDTH=$DMA_WIDTH
solution file set ../inc/ad03_cxx_catapult.hpp -args -DDMA_WIDTH=$DMA_WIDTH
solution file set ../tb/main.cpp -args -DDMA_WIDTH=$DMA_WIDTH

#
# Output
#

# Verilog only
solution option set Output/OutputVHDL false
solution option set Output/OutputVerilog true

# Package output in Solution dir
solution option set Output/PackageOutput true
solution option set Output/PackageStaticFiles true

# Add Prefix to library and generated sub-blocks
solution option set Output/PrefixStaticFiles true
solution options set Output/SubBlockNamePrefix "esp_acc_${ACCELERATOR}_"

# Do not modify names
solution option set Output/DoNotModifyNames true

go new

#
#
#

go analyze

#
#
#


# Set the top module and inline all of the other functions.

# 10.4c
#directive set -DESIGN_HIERARCHY ${ACCELERATOR}

# 10.5
solution design set $ACCELERATOR -top

#directive set PRESERVE_STRUCTS false

#
#
#

go compile

# Run C simulation.
if {$opt(csim)} {
    flow run /SCVerify/launch_make ./scverify/Verify_orig_cxx_osci.mk {} SIMTOOL=osci sim
}

#
#
#

# Run HLS.
if {$opt(hsynth)} {

    if {[lsearch $fpga_techs $TECH] >= 0} {
	solution library \
	    add mgc_Xilinx-$FPGA_FAMILY$FPGA_SPEED_GRADE\_beh -- \
	    -rtlsyntool Vivado \
	    -manufacturer Xilinx \
	    -family $FPGA_FAMILY \
	    -speed $FPGA_SPEED_GRADE \
	    -part $FPGA_PART_NUM
	solution library add Xilinx_RAMS
	solution library add Xilinx_ROMS
	solution library add Xilinx_FIFO
    } else {
	solution library add nangate-45nm_beh -- -rtlsyntool DesignCompiler -vendor Nangate -technology 045nm
    }

    # For Catapult 10.5: disable all sequential clock-gating
    directive set GATE_REGISTERS false

    go libraries

    #
    #
    #

    directive set -CLOCKS { \
        clk { \
            -CLOCK_PERIOD 4.0 \
            -CLOCK_EDGE rising \
            -CLOCK_HIGH_TIME 2.0 \
            -CLOCK_OFFSET 0.000000 \
            -CLOCK_UNCERTAINTY 0.0 \
            -RESET_KIND sync \
            -RESET_SYNC_NAME rst \
            -RESET_SYNC_ACTIVE low \
            -RESET_ASYNC_NAME arst_n \
            -RESET_ASYNC_ACTIVE low \
            -ENABLE_NAME {} \
            -ENABLE_ACTIVE high \
        } \
    }

    # BUGFIX: This prevents the creation of the empty module CGHpart. In the
    # next releases of Catapult HLS, this may be fixed.
    directive set /$ACCELERATOR -GATE_EFFORT normal

    go assembly

    #
    #
    #

    # Top-Module I/O
    directive set /$ACCELERATOR/conf_info:rsc -MAP_TO_MODULE ccs_ioport.ccs_in_wait
    directive set /$ACCELERATOR/dma_read_ctrl:rsc -MAP_TO_MODULE ccs_ioport.ccs_out_wait
    directive set /$ACCELERATOR/dma_write_ctrl:rsc -MAP_TO_MODULE ccs_ioport.ccs_out_wait
    directive set /$ACCELERATOR/dma_read_chnl:rsc -MAP_TO_MODULE ccs_ioport.ccs_in_wait
    directive set /$ACCELERATOR/dma_write_chnl:rsc -MAP_TO_MODULE ccs_ioport.ccs_out_wait
    directive set /$ACCELERATOR/acc_done:rsc -MAP_TO_MODULE ccs_ioport.ccs_sync_out_vld

    # Arrays
    #directive set /$ACCELERATOR/core/plm_in.data:rsc -MAP_TO_MODULE Xilinx_RAMS.BLOCK_1R1W_RBW
    #directive set /$ACCELERATOR/core/plm_out.data:rsc -MAP_TO_MODULE Xilinx_RAMS.BLOCK_1R1W_RBW

    directive set /$ACCELERATOR/core/plm_in.data:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/plm_out.data:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/compute_wrapper<plm_in_t,plm_out_t>:tmp:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/anomaly_detector:layer2_out:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/anomaly_detector:layer4_out:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/anomaly_detector:layer5_out:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/anomaly_detector:layer6_out:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/anomaly_detector:layer8_out:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/anomaly_detector:layer9_out:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/anomaly_detector:layer10_out:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/anomaly_detector:layer12_out:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/anomaly_detector:layer13_out:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/anomaly_detector:layer14_out:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/anomaly_detector:layer16_out:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/anomaly_detector:layer17_out:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/anomaly_detector:layer18_out:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/anomaly_detector:layer20_out:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/anomaly_detector:layer21_out:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/nnet::dense_resource_rf_gt_nin_rem0<input_t,layer2_t,config2>:acc:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/nnet::dense_resource_rf_gt_nin_rem0<layer5_t,layer6_t,config6>:acc:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/nnet::dense_resource_rf_gt_nin_rem0<layer9_t,layer10_t,config10>:acc:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/nnet::dense_resource_rf_gt_nin_rem0<layer13_t,layer14_t,config14>:acc:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/nnet::dense_resource_rf_gt_nin_rem0<layer17_t,layer18_t,config18>:acc:rsc -MAP_TO_MODULE {[Register]}
    directive set /$ACCELERATOR/core/nnet::dense_resource_rf_gt_nin_rem0<layer21_t,layer22_t,config22>:acc:rsc -MAP_TO_MODULE {[Register]}

    directive set /ad03_cxx_catapult/w2.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/b4.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/w14.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/w22.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/s4.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/w6.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/s8.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/b8.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/w10.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/s12.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/b12.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/s16.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/b16.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/w18.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/s20.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/b20.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/nnet::dense_resource_rf_gt_nin_rem0<layer9_t,layer10_t,config10>:acc.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/nnet::dense_resource_rf_gt_nin_rem0<layer17_t,layer18_t,config18>:acc.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/nnet::dense_resource_rf_gt_nin_rem0<layer21_t,layer22_t,config22>:acc.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/nnet::dense_resource_rf_gt_nin_rem0<input_t,layer2_t,config2>:outidx.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/nnet::dense_resource_rf_gt_nin_rem0<layer5_t,layer6_t,config6>:outidx.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/nnet::dense_resource_rf_gt_nin_rem0<layer9_t,layer10_t,config10>:outidx.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/nnet::dense_resource_rf_gt_nin_rem0<layer13_t,layer14_t,config14>:outidx.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/nnet::dense_resource_rf_gt_nin_rem0<layer17_t,layer18_t,config18>:outidx.rom:rsc -MAP_TO_MODULE {[Register]}
    directive set /ad03_cxx_catapult/nnet::dense_resource_rf_gt_nin_rem0<layer21_t,layer22_t,config22>:outidx.rom:rsc -MAP_TO_MODULE {[Register]}

    # Loops
    # 1 function
    #directive set /$ACCELERATOR/core/BATCH_LOOP -PIPELINE_INIT_INTERVAL 1
    #directive set /$ACCELERATOR/core/BATCH_LOOP -PIPELINE_STALL_MODE flush

    # Loops performance tracing


    # Area vs Latency Goals
    directive set /$ACCELERATOR/core -EFFORT_LEVEL high
    directive set /$ACCELERATOR/core -DESIGN_GOAL Latency

    if {$opt(debug) != 1} {
        go architect

        #
        #
        #

        go allocate

        #
        # RTL
        #

        #directive set ENABLE_PHYSICAL true

        go extract

        #
        #
        #

        if {$opt(rtlsim)} {
            #flow run /SCVerify/launch_make ./scverify/Verify_concat_sim_${ACCELERATOR}_v_msim.mk {} SIMTOOL=msim sim
            flow run /SCVerify/launch_make ./scverify/Verify_concat_sim_${ACCELERATOR}_v_msim.mk {} SIMTOOL=msim simgui
        }

        if {$opt(lsynth)} {
            flow run /Vivado/synthesize -shell vivado_concat_v/concat_${ACCELERATOR}.v.xv
            #flow run /Vivado/synthesize vivado_concat_v/concat_${ACCELERATOR}.v.xv
        }
    }
}

project save

puts "***************************************************************"
puts "uArch: Single block"
puts "***************************************************************"
puts "Done!"

