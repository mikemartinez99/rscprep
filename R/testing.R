test <- "/Volumes/rc/lab/G/GSR_Active/Labs/Pioli/260625_10X_Flex/analysis/"

manifest <- data.frame(
    Sample = c("Pool_1"),
    Path = "/Volumes/rc/lab/G/GSR_Active/Labs/Pioli/260625_10X_Flex/analysis/Pool1/outs/",
    multiplexed = TRUE
)

#----- Testing function with both manifest and path
seuratFrom10x(manifest = manifest)
seuratFrom10x(path = "/Volumes/rc/lab/G/GSR_Active/Labs/Pioli/260625_10X_Flex/analysis/Pool1/outs/", multiplexed = TRUE)
