echo "run lasso"
Rscript run_rf_model_notime.R
wait
echo "run random forest"
Rscript run_lasso_model_notime.R
wait
echo "I am using the lasso features to test different models"
Rscript run_compare_multiple_models_v2.R
wait
echo "I am using the RF features to test different models"
Rscript run_explore_res_rf_model_notime.R
wait
echo "I am comparing lasso with RF"
Rscript run_compare_lasso_with_RF.R
wait
echo "I am comparing lasso with RF2"
Rscript run_models_common_genes_lasso_rf.R
