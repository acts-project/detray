from detray_sympy.common import (
    cxx_printer,
    my_expression_print,
)


def gen_cxx_code(function_name, inputs, outputs, run_cse=True):
    printer = cxx_printer

    template_types = []

    for i, j in inputs:
        template_types.append("%s_t" % i)
    for i, j in outputs:
        template_types.append("%s_t" % i)

    lines = []

    lines.append(
        "template <%s>" % (", ".join("typename %s" % s for s in template_types))
    )
    lines.append("DETRAY_HOST_DEVICE void inline %s (" % function_name)
    lines.append(", ".join("const %s_t & %s" % (i, i) for i, _ in inputs) + ",")
    lines.append(", ".join("%s_t & %s" % (i, i) for i, _ in outputs))
    lines.append(") {")

    code = my_expression_print(
        printer,
        outputs,
        [x[0] for x in outputs],
        run_cse=run_cse,
    )
    lines.extend([f"  {l}" for l in code.split("\n")])

    lines.append("}")

    return "\n".join(lines)
