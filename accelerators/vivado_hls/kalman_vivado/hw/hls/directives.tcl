# Copyright (c) 2011-2023 Columbia University, System Level Design Group
# SPDX-License-Identifier: Apache-2.0

# User-defined configuration ports
# <<--directives-param-->>
set_directive_interface -mode ap_none "top" conf_info_iter
set_directive_interface -mode ap_none "top" conf_info_x_dim
set_directive_interface -mode ap_none "top" conf_info_z_dim

# Insert here any custom directive
#set_directive_dataflow "top/go"
set_directive_array_partition -type cyclic -factor 2 -dim 1 "top" _inbuff
set_directive_array_partition -type cyclic -factor 2 -dim 1 "top" _outbuff

set_directive_array_partition -type cyclic -factor 2 -dim 1 "load" tmp
set_directive_pipeline "load/load_label1"

# set_directive_pipeline "compute"
# set_directive_unroll -off "compute/LOOP_SR_1"
