using DelimitedFiles
using Plots

# Read Ramsey data
# Load Ramsey 01 data
if (data_set_option=="short")
    ramsey_01_0_file = "data-set-20220413/population_ramsey_01_half_period_1000000.0_1000_20220413_4000_5_0.txt"
    ramsey_01_1_file = "data-set-20220413/population_ramsey_01_half_period_1000000.0_1000_20220413_4000_5_1.txt"
    ramsey_01_2_file = "data-set-20220413/population_ramsey_01_half_period_1000000.0_1000_20220413_4000_5_2.txt"
    t_ramsey_01_file = "data-set-20220413/darktime_ramsey_01_half_period_1000000.0_1000_20220413_4000_5.txt"

    ramsey_12_0_file = "data-set-20220413/population_ramsey_12_1000000.0_4000_5_1000_confusion_0.txt"
    ramsey_12_1_file = "data-set-20220413/population_ramsey_12_1000000.0_4000_5_1000_confusion_1.txt"
    ramsey_12_2_file = "data-set-20220413/population_ramsey_12_1000000.0_4000_5_1000_confusion_2.txt"
    t_ramsey_12_file = "data-set-20220413/darktime_ramsey_12_1000000.0_5_4000_1000_20220413_confusion.txt"
else
    ramsey_01_0_file = "data-set-20220413/population_ramsey_01_half_period_1000000.0_1000_20220413_80000_20_0.txt"
    ramsey_01_1_file = "data-set-20220413/population_ramsey_01_half_period_1000000.0_1000_20220413_80000_20_1.txt"
    ramsey_01_2_file = "data-set-20220413/population_ramsey_01_half_period_1000000.0_1000_20220413_80000_20_2.txt"
    t_ramsey_01_file = "data-set-20220413/darktime_ramsey_01_half_period_1000000.0_1000_20220413_80000_20.txt"

    ramsey_12_0_file = "data-set-20220413/population_ramsey_12_1000000.0_80000_20_1000_confusion_0.txt"
    ramsey_12_1_file = "data-set-20220413/population_ramsey_12_1000000.0_80000_20_1000_confusion_1.txt"
    ramsey_12_2_file = "data-set-20220413/population_ramsey_12_1000000.0_80000_20_1000_confusion_2.txt"
    t_ramsey_12_file = "data-set-20220413/darktime_ramsey_12_1000000.0_20_80000_1000_20220413_confusion.txt"
end