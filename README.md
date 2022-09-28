# ZF_expression_browser

> Zebrafish Gene Expression Browser (eggs, development)

***How to use this Shiny App:***

1. Create a directory, e.g. ZFApp
   - Put the files app.R and allAppData.R in this directory
2. Install R
3. Install shiny and DT libraries. In R type:
   - install.packages('shiny')
   - install.packages('DT')
4. load the shiny library into your R session:
   - library(shiny)
5. Issue the runApp command from within R while pointing at the directory you created in 1:
   - runApp('ZFApp')

***Alternatively, if you have git installed:***

1. Install R
2. Install shiny and DT libraries. In R type:
   - install.packages('shiny')
   - install.packages('DT')
3. Open your git bash shell
   -  git clone https://github.com/rauwerda/ZF_expression_browser.git
4. load the shiny library into your R session:
   - library(shiny)
5. Issue the runApp command from within R while pointing at the directory you created in 1:
   - runApp("ZF_expression_browser")
