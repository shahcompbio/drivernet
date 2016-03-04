## Define accessors for the DriverNetResult class object.
setMethod("drivers", "DriverNetResult", function(x) x@drivers)

setMethod("actualEvents", "DriverNetResult", function(x) x@actualEvents)

setMethod("totalEvents", "DriverNetResult", function(x) x@totalEvents)
