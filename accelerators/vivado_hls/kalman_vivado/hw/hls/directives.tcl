# Copyright (c) 2011-2023 Columbia University, System Level Design Group
# SPDX-License-Identifier: Apache-2.0

# User-defined configuration ports
# <<--directives-param-->>
set_directive_interface -mode ap_none "top" conf_info_iter
set_directive_interface -mode ap_none "top" conf_info_x_dim
set_directive_interface -mode ap_none "top" conf_info_z_dim

# Insert here any custom directive
#set_directive_dataflow "top/go"
#set_directive_array_partition -type cyclic -factor 2 -dim 1 "top" _inbuff
set_directive_array_partition -type cyclic -factor 2 -dim 1 "top" _outbuff

set_directive_array_partition -type complete -dim 0 "load" tmp
set_directive_pipeline "load/load_label1"

# set_directive_pipeline "compute"
# set_directive_unroll -off "compute/LOOP_SR_1"

# set_directive_dataflow "qr_inverse_top"
set_directive_inline -recursive "qr_inverse_top"
#set_directive_inline -recursive "matrix_multiply_top"

set_directive_pipeline "load/load_pred"
# set_directive_pipeline "compute/LOOP_INT2_1"
# set_directive_array_partition -type complete -dim 2 "compute" inter2
# set_directive_array_partition -type complete -dim 0 "top" X
# set_directive_array_partition -type complete -dim 0 "top" P
# set_directive_array_map -instance ALL -mode horizontal compute X
# set_directive_array_map -instance ALL -mode horizontal compute P
# set_directive_array_map -instance ALL -mode horizontal compute H
# set_directive_array_map -instance ALL -mode horizontal compute R
# set_directive_array_map -instance ALL -mode horizontal compute Z
# set_directive_array_map -instance ALL -mode horizontal compute F
# set_directive_array_map -instance ALL -mode horizontal compute Q
# set_directive_array_map -instance ALL -mode horizontal compute X_pred
# set_directive_array_map -instance ALL -mode horizontal compute P_pred

