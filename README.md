# Metrics-Project

We estimate the static effect on the outcomes of interest of changes in inputs tariff
between the years 1996 and 1999, as such period is characterised by district which
experience a change in input tariff (switchers) and district which do not (stayers) which is
the fundamental assumption of the diff in diff procedure.

did_multiplegt_dyn(
    df = indonesia,
    outcome = "p0",
    group = "distid",
    time = "year",
    effects=1,
    treatment = "TarIn",
    normalized = T
)
