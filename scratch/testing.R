# Scratch: manual checks against real data on the lab volume.
# Lives outside R/ on purpose -- anything in R/ is sourced by devtools::load_all(),
# which would make every package load (and every test run) depend on this share
# being mounted. Run it the usual way: devtools::load_all() then these lines.

devtools::load_all()

test <- "/Volumes/rc/lab/G/GSR_Active/Labs/Pioli/260625_10X_Flex/analysis/"

manifest <- data.frame(
    Sample = c("Pool_1"),
    Path = "/Volumes/rc/lab/G/GSR_Active/Labs/Pioli/260625_10X_Flex/analysis/Pool1/outs/",
    multiplexed = TRUE
)

#----- Testing function with both manifest and path
seuratFrom10x(manifest = manifest)
seuratFrom10x(path = "/Volumes/rc/lab/G/GSR_Active/Labs/Pioli/260625_10X_Flex/analysis/Pool1/outs/", multiplexed = TRUE)
