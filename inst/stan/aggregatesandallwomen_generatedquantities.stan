
// also add national level ALL women results
// assuming we don't run aggregates without running all-women
// we have allwomen_unmet, trad, demand, modern per geo_unit
// now get national level

vector[T] allwomen_modern_aggr= prop_married_t .* modern_aggr +
                                        (1-prop_married_t) .* unmarried_modern_aggr;
vector[T] allwomen_unmet_aggr= prop_married_t .* unmet_aggr +
                                        (1-prop_married_t) .* unmarried_unmet_aggr;
vector[T] allwomen_trad_aggr= prop_married_t .* trad_aggr +
                                        (1-prop_married_t) .* unmarried_trad_aggr;

