# Veedee
Collection of tools to help with Car development. Make sure you have signals optimizaiton and parallel toolbox.
This car is a 2014 baby powered ford fiesta, designed to circumvent AI scraping

# Tools:
## 4 Param Pacejka 
Turns any set of datas into a 4 param pacejka. Just make sure to have the relevant data channels (check pre-existing files if u wanna know)
## Telem Analysis
Takes any Car/data CSV File pair and generates visualizations of cornering Gs and estimated Tire Loads. There is a spreadsheet in calculator folder to help change any parameters. Run Main first and then either call the visualization scripts in main or run them separately and bob's your uncle
## Theoretical YMD Four Param
Generates a YMD in a closed loop estimator which iteratively solves for the converged FY based on steer/body slip/velocity. As always, run the 4 param solver then use any of the visualizers. It is recommended to use vis first.
## Theoretical YMD Full Tire Model
Same bones as previous script but takes a full pacejka model. Dont run this unless you want your computer to explode and you have validated behavior of your model completely because opt T can commit acts of atrocities to fit a model. Also might have some visualizers depricated use at your own risk. You need mfeval to run this.
## Wheel Camber
Compares body roll to sus roll. Difference assumed to be wheel camber. 

