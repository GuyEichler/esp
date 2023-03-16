# Copyright (c) 2011-2022 Columbia University, System Level Design Group
# SPDX-License-Identifier: Apache-2.0

# User-defined configuration ports
# <<--directives-param-->>
set_directive_interface -mode ap_none "top" conf_info_avg
set_directive_interface -mode ap_none "top" conf_info_key_length
set_directive_interface -mode ap_none "top" conf_info_std
set_directive_interface -mode ap_none "top" conf_info_R
set_directive_interface -mode ap_none "top" conf_info_L
set_directive_interface -mode ap_none "top" conf_info_key_batch
set_directive_interface -mode ap_none "top" conf_info_key_num
set_directive_interface -mode ap_none "top" conf_info_val_num

# Insert here any custom directive
set_directive_loop_tripcount -min 256 -max 256 -avg 256 "top/go_2"
#set_directive_unroll -factor 2 "store_val/store_label1"
# set_directive_dataflow "top/go_2"

set_directive_pipeline "compute_val/ASSIGN_LOOP"
set_directive_pipeline -II 3 "compute/COMPUTE_LOOP"
