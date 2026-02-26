import sys
import sympy

from detray_sympy.common import name_expr
from detray_sympy.output import write_out_file
from detray_sympy.codegen import gen_cxx_code


def gen_code():
    transport_jacobian = (
        sympy.MatrixSymbol("transport_jacobian", 8, 8).as_explicit().as_mutable()
    )

    b2f_dpos_dloc = sympy.MatrixSymbol("b2f_dpos_dloc", 3, 2).as_explicit().as_mutable()
    b2f_ddir_dangle = (
        sympy.MatrixSymbol("b2f_ddir_dangle", 3, 2).as_explicit().as_mutable()
    )
    b2f_ddir_dangle[2, 0] = 0
    b2f_dpos_dangle = (
        sympy.MatrixSymbol("b2f_dpos_dangle", 3, 2).as_explicit().as_mutable()
    )

    path_to_free_derivative = (
        sympy.MatrixSymbol("path_to_free_derivative", 8, 1).as_explicit().as_mutable()
    )
    path_to_free_derivative[3, 0] = 0

    free_to_path_derivative = (
        sympy.MatrixSymbol("free_to_path_derivative", 1, 8).as_explicit().as_mutable()
    )
    free_to_path_derivative[0, 3] = 0
    free_to_path_derivative[0, 7] = 0

    f2b_dloc_dpos = sympy.MatrixSymbol("f2b_dloc_dpos", 2, 3).as_explicit().as_mutable()
    f2b_dangle_ddir = (
        sympy.MatrixSymbol("f2b_dangle_ddir", 2, 3).as_explicit().as_mutable()
    )
    f2b_dangle_ddir[0, 2] = 0

    bound_to_free_jacobian = sympy.zeros(8, 6)
    bound_to_free_jacobian[0:3, 0:2] = b2f_dpos_dloc
    bound_to_free_jacobian[0:3, 2:4] = b2f_dpos_dangle
    bound_to_free_jacobian[4:7, 2:4] = b2f_ddir_dangle
    bound_to_free_jacobian[7, 4] = 1
    bound_to_free_jacobian[3, 5] = 1

    free_to_bound_jacobian = sympy.zeros(6, 8)
    free_to_bound_jacobian[0:2, 0:3] = f2b_dloc_dpos
    free_to_bound_jacobian[2:4, 4:7] = f2b_dangle_ddir
    free_to_bound_jacobian[4, 7] = 1
    free_to_bound_jacobian[5, 3] = 1

    full_jacobian = (
        free_to_bound_jacobian
        * (path_to_free_derivative * free_to_path_derivative + sympy.eye(8))
        * transport_jacobian
        * bound_to_free_jacobian
    )

    tmp = sympy.eye(6)
    tmp[0:4, 0:5] = full_jacobian[0:4, 0:5]
    tmp[5:6, 0:5] = full_jacobian[5:6, 0:5]
    full_jacobian = tmp

    input_name_exprs = [
        name_expr("transport_jacobian", transport_jacobian),
        name_expr("b2f_dpos_dloc", b2f_dpos_dloc),
        name_expr("b2f_ddir_dangle", b2f_ddir_dangle),
        name_expr("b2f_dpos_dangle", b2f_dpos_dangle),
        name_expr("path_to_free_derivative", path_to_free_derivative),
        name_expr("free_to_path_derivative", free_to_path_derivative),
        name_expr("f2b_dloc_dpos", f2b_dloc_dpos),
        name_expr("f2b_dangle_ddir", f2b_dangle_ddir),
    ]
    output_name_exprs = [name_expr("full_jacobian", full_jacobian)]
    code = gen_cxx_code(
        "update_full_jacobian_impl",
        input_name_exprs,
        output_name_exprs,
        run_cse=True,
    )
    return code


if __name__ == "__main__":
    if len(sys.argv) > 1:
        output = sys.argv[1]
    else:
        output = None

    code = gen_code()

    write_out_file(code, output)
