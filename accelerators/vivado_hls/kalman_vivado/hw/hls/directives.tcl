# Copyright (c) 2011-2023 Columbia University, System Level Design Group
# SPDX-License-Identifier: Apache-2.0

# User-defined configuration ports
# <<--directives-param-->>
set_directive_interface -mode ap_none "top" conf_info_iter
set_directive_interface -mode ap_none "top" conf_info_x_dim
set_directive_interface -mode ap_none "top" conf_info_z_dim

# Insert here any custom directive
# set_directive_dataflow "top/go"
#set_directive_dataflow "compute"
#set_directive_array_partition -type cyclic -factor 2 -dim 1 "top" _inbuff
set_directive_array_partition -type cyclic -factor 2 -dim 1 "top" _outbuff

set_directive_array_partition -type complete -dim 0 "load" tmp
# set_directive_pipeline "load/load_label1"
set_directive_pipeline "load/load_init1"
set_directive_dependence -variable X -type inter -dependent false "load/load_init1"
set_directive_dependence -variable X_pred -type inter -dependent false "load/load_init1"
set_directive_dependence -variable Z -type inter -dependent false "load/load_init1"
set_directive_dependence -variable P -type inter -dependent false "load/load_init1"
set_directive_dependence -variable P_pred -type inter -dependent false "load/load_init1"
set_directive_dependence -variable F -type inter -dependent false "load/load_init1"
set_directive_dependence -variable Q_kal -type inter -dependent false "load/load_init1"
set_directive_dependence -variable H -type inter -dependent false "load/load_init1"
set_directive_dependence -variable R_kal -type inter -dependent false "load/load_init1"


# set_directive_pipeline "compute"
# set_directive_unroll -off "compute/LOOP_SR_1"

# set_directive_dataflow "qr_inverse_top"
# set_directive_inline -recursive "qr_inverse_top"
#set_directive_inline -recursive "matrix_multiply_top"

set_directive_pipeline "load/load_pred"
set_directive_dependence -variable X -type inter -dependent false "load/load_pred"
set_directive_dependence -variable P -type inter -dependent false "load/load_pred"
# set_directive_pipeline "compute/LOOP_INT2_2"
set_directive_pipeline "compute/LOOP_S_inv_2"
set_directive_dependence -variable S_inv -type inter -dependent false "compute/LOOP_S_inv_2"
set_directive_pipeline -II 1 "compute/LOOP_SR_2"
set_directive_dependence -variable S -type inter -dependent false "compute/LOOP_SR_2"
set_directive_pipeline "compute/LOOP_Y"
set_directive_dependence -variable Y -type inter -dependent false "compute/LOOP_Y"
set_directive_pipeline "compute/LOOP_X_PRED"
set_directive_dependence -variable X_pred -type inter -dependent false "compute/LOOP_X_PRED"
set_directive_pipeline "compute/LOOP_INT7_2"
set_directive_dependence -variable inter7 -type inter -dependent false "compute/LOOP_INT7_2"
set_directive_pipeline "compute/LOOP_OUT"
set_directive_dependence -variable _outbuff -type inter -dependent false "compute/LOOP_OUT"
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

# set_directive_array_partition -type cyclic -factor 2 -dim 2 compute S
# set_directive_array_partition -type cyclic -factor 2 -dim 2 compute S_inv
# set_directive_array_partition -type block -factor 4 -dim 1 compute S_inv
# set_directive_array_map -instance RS -mode horizontal compute R
# set_directive_array_map -instance RS -mode horizontal compute S_inv

# set_directive_array_map -instance inter13 -mode horizontal compute inter3
# set_directive_array_map -instance inter13 -mode horizontal compute inter1
# set_directive_array_map -instance inter26 -mode horizontal compute inter2
# set_directive_array_map -instance inter26 -mode horizontal compute inter6

# set_directive_array_map -instance inter1346 -mode horizontal compute inter16
# set_directive_array_map -instance inter1346 -mode horizontal compute inter34
