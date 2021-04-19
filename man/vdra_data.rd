\name{vdra_data}
\alias{vdra_data}
\title{
    Simulated data from a weight loss study.
}
\description{
    Simulated data based on a weight loss study performed by Harvard School of Medicine.  The original study was looking for comorbidities to predict weight loss / weight gain in obese subjects.  The various outcomes in this simulated data set don't really correlate to each other.  The purpose of this simulated data set is to illustrate the funcationality of the package, not to draw valid statistical inferences.
}
\usage{
    vdra_data
}
\format{
  \tabular{ll}{
    Change_BMI:\tab Continuous response used for linear regression.\cr
    WtLost:\tab Binary response used for logistic regression.\cr
    Time:\tab Used for time to event in Cox regression.  Continuous.  Ranges from 1 to 459.\cr
    Status:\tab Used for censoring in Cox regression.  A binary categroical variable.\cr
    Exposure:\tab A binary categorical variable.\cr
    Age:\tab A continuous variable ranging from 3 to 80.\cr
    ComorbidScore:\tab A ordinal variable with 11 variables: 0 to 10.\cr
    NumRx:\tab A ordinal varible with 15 levels:  -2 to 12.\cr
    BMI_pre: \tab A continuous variable ranging from 35.01 to 92.79.\cr
    Race:\tab A factor with 6 levels: "Race 0" to "Race 5". \cr
    Sex:\tab A binary factor with 2 levels: M and F.\cr
  }
}
