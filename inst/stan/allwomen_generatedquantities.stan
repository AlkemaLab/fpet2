// allwomen results
// we have unmet, trad, and d_Eta (demand), for "" and "unmarried_"
// results per geounit
matrix[n_geounit, T] allwomen_unmet = prop_married_rt .* unmet +
                                        (1-prop_married_rt) .* unmarried_unmet;
matrix[n_geounit, T] allwomen_trad = prop_married_rt .* trad +
                                        (1-prop_married_rt) .* unmarried_trad;
matrix[n_geounit, T] allwomen_demand = prop_married_rt .* d_Eta +
                                        (1-prop_married_rt) .* unmarried_d_Eta;
matrix[n_geounit, T] allwomen_modern = allwomen_demand - allwomen_unmet;



