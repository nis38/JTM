## Boyce Index ##

hsm10 <- hsmmod(clim2010, occur00s)[[1]]
plot(hsm10)
boyc10 <- ecospat.boyce(hsm10, occur00s[c(1:75, 77, 78),], method = "spearman")
boyc10$cor

hsm11 <- hsmmod(clim2011, occur00s)[[1]]
plot(hsm11)
boyc11 <- ecospat.boyce(hsm11, occur00s[c(1:75, 77, 78),], method = "spearman")
boyc11$cor

hsm12 <- hsmmod(clim2012, occur00s)[[1]]
plot(hsm12)
boyc12 <- ecospat.boyce(hsm12, occur00s[c(1:75, 77, 78),], method = "spearman")
boyc12$cor

hsm13 <- hsmmod(clim2013, occur00s)[[1]]
plot(hsm13)
boyc13 <- ecospat.boyce(hsm13, occur00s[c(1:75, 77, 78),], method = "spearman")
boyc13$cor

hsm14 <- hsmmod(clim2014, occur00s)[[1]]
plot(hsm14)
boyc14 <- ecospat.boyce(hsm14, occur00s[c(1:75, 77, 78),], method = "spearman")
boyc14$cor

hsm15 <- hsmmod(clim2015, occur00s)[[1]]
plot(hsm15)
boyc15 <- ecospat.boyce(hsm15, occur00s[c(1:75, 77, 78),], method = "spearman")
boyc15$cor

raster::extract(hsm16, occur00s)

boyc16 <- ecospat.boyce(hsm16, occur00s[c(1:75, 77, 78),], method = "spearman")
boyc16$cor


raster::extract(hsm17, occur00s)

boyc17 <- ecospat.boyce(hsm17, occur00s[c(1:75, 77, 78),], method = "spearman")
boyc17$cor



boyc18 <- ecospat.boyce(hsm18, occur00s[c(1:75, 77, 78),], method = "spearman")
boyc18$cor


mean(boyc10$cor, boyc11$cor, boyc12$cor, boyc13$cor, boyc14$cor, boyc15$cor,
     boyc16$cor, boyc17$cor, boyc18$cor)
sd(c(boyc10$cor, boyc11$cor, boyc12$cor, boyc13$cor, boyc14$cor, boyc15$cor,
    boyc16$cor, boyc17$cor, boyc18$cor))
