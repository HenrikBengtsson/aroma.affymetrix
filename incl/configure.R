invisible({
# buildScripts/
file.rename("inst/bu", "inst/buildScripts")

# testScripts/
file.rename("inst/te/ad", "inst/t/addons")
file.rename("inst/te/ar/sy/ch", "inst/t/ar/sy/chipTypes")
file.rename("inst/te/ar/sy", "inst/t/ar/system")
file.rename("inst/te/ar", "inst/t/archive")
file.rename("inst/te/co", "inst/t/complete")
file.rename("inst/te/sy/ch", "inst/t/system/chipTypes")
file.rename("inst/te/sy", "inst/t/system")
file.rename("inst/te", "inst/testScripts")
})

