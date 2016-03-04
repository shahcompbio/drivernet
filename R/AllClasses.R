## Class to store result returned by "computeDrivers" function
setClass("DriverNetResult",
         representation=representation(
           drivers="character",
           actualEvents="list",
           totalEvents="numeric"))