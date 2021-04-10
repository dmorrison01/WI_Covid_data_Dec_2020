# WI_Covid_data_Dec_2020
repository for extension of Hybrid Shewhart chart functions

I generalized the names of variables used in the Hybrid Shewhart chart files.  In the previous work [IHI_Covid_Public](https://github.com/klittle314/IHI_Covid_Public) many of the names of variables in data frames explicitly referred to DEATHS.  Otherwise, the logic is the same as described in the IHI_Covid_Public [README.md file](https://github.com/klittle314/IHI_Covid_Public/blob/main/README.md).  Examples use data from the Wisconsin Department of Public Health.

The core function find_phase_dates is generalized to include a parameter Epoch3_4_transition, which allows the user to set a threshold in the transition from a phase in Epoch 3 to a phase in Epoch 4.  For death series, our original target, we set the threshold as 2 deaths, which translated into a requirement that the calculated lower limit be less than 2 deaths to consider a transition from Epoch 3 to Epoch 4.

I include an RMarkdown file and output HTML format: [markdown](https://github.com/klittle314/WI_Covid_data_Dec_2020/blob/master/Wisc%20data%20Dec%202020.Rmd) to generate interactive plots (set output_format parameter to HTML and knit to HTML) via plotly.

I also include a summary of World Covid attributed deaths in a markdown file, as another illustration [world markdown](https://github.com/klittle314/WI_Covid_data_Dec_2020/blob/master/world_output.Rmd)
