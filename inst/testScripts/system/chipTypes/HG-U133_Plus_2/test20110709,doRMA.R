library("aroma.affymetrix");

ces <- doRMA("Affymetrix-HeartBrain", chipType="HG-U133_Plus_2", verbose=-8);
print(ces);
