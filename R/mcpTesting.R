#### Testing regression with multiple change points ####

### Complete Skull
PCA <- Croc_Lateral_PCA
CACs <- Croc_Lateral_CACs

LM_MCP_array <- data.frame(PCA$x,CACs, log(Croc_Lateral_GPA$Csize), lateral_covariates$species, lateral_covariates$Ecomorph, lateral_covariates$simple_ecomorph, stringsAsFactors=FALSE)
dimnames(LM_MCP_array)[[2]][64] <- "logCS"
dimnames(LM_MCP_array)[[2]][65] <- "species"
dimnames(LM_MCP_array)[[2]][66] <- "ecomorph"
dimnames(LM_MCP_array)[[2]][67] <- "simple_ecomorph"

#Null model - no break points and single slope/intercept
Null_CAC_model = list(
  CACs ~ logCS
)

Null_CAC_LM_mcp <- mcp(Null_CAC_model, LM_MCP_array)
plot(Null_CAC_LM_mcp)
summary(Null_CAC_LM_mcp)
#

#model for single break point with multiple slopes
CAC_LM_model = list(
  CACs ~ 1 + logCS,
  ~ 0 + logCS
)

CAC_LM_mcp <- mcp(CAC_LM_model, LM_MCP_array)
plot(CAC_LM_mcp)
summary(CAC_LM_mcp)
#

#model for single break point with multiple intercepts and slopes
CAC_LM_int_model = list(
  CACs ~ 1 + logCS,
  ~ 1 + logCS
)

CAC_LM_int_mcp <- mcp(CAC_LM_int_model, LM_MCP_array)
plot(CAC_LM_int_mcp)
summary(CAC_LM_int_mcp)
#

#model for single break point with variance among species and multiple intercepts and slopes in both lines
CAC_LM_sp_model_V2 = list(
  CACs ~ 1 + (1|species) + logCS,
  ~ (0|species) + logCS
)

CAC_LM_sp_mcp_V2 <- mcp(CAC_LM_sp_model_V2, LM_MCP_array)
plot(CAC_LM_sp_mcp_V2)
plot(CAC_LM_sp_mcp_V2, facet_by = "species")
summary(CAC_LM_sp_mcp_V2)
#

#model for single break point with variance among species with multiple intercepts and slopes for only the second line
CAC_LM_sp_model = list(
  CACs ~ 1 + logCS,
  1 + (1|species) ~ (0|species) + logCS
)

CAC_LM_sp_mcp <- mcp(CAC_LM_sp_model, LM_MCP_array)
plot(CAC_LM_sp_mcp)
plot(CAC_LM_sp_mcp, facet_by = "species")
summary(CAC_LM_sp_mcp)
ranef(CAC_LM_sp_mcp)
#

#model for single break point with variance among species and multiple intercepts and slopes in both lines
CAC_LM_eco_model_V2 = list(
  CACs ~ 1 + (1|simple_ecomorph) + logCS,
  ~ (0|simple_ecomorph) + logCS
)

CAC_LM_eco_mcp_V2 <- mcp(CAC_LM_eco_model_V2, LM_MCP_array)
plot(CAC_LM_eco_mcp_V2)
summary(CAC_LM_eco_mcp_V2)


#model for single break point with variance among species and multiple intercepts and slopes in both lines
CAC_LM_eco_model_V3 = list(
  CACs ~ (1|simple_ecomorph) + logCS,
  1 + (1|simple_ecomorph) ~ (0|simple_ecomorph) + logCS
)

CAC_LM_eco_mcp_V3 <- mcp(CAC_LM_eco_model_V3, LM_MCP_array)
plot(CAC_LM_eco_mcp_V3)
plot(CAC_LM_eco_mcp_V3, facet_by = "simple_ecomorph")
summary(CAC_LM_eco_mcp_V3)

#model for single break point with variance among species and multiple intercepts and slopes in the second line only
CAC_LM_eco_model = list(
  CACs ~ 1 + logCS,
  1 + (1|simple_ecomorph) ~ (0|simple_ecomorph) + logCS
)

CAC_LM_eco_mcp <- mcp(CAC_LM_eco_model, LM_MCP_array)
plot(CAC_LM_eco_mcp)
plot(CAC_LM_eco_mcp, facet_by = "simple_ecomorph")
summary(CAC_LM_eco_mcp)
ranef(CAC_LM_eco_mcp)
#

#model for two break points with multiple intercepts and slopes
CAC_2breaks_int_model = list(
  CACs ~ 1 + logCS,
  ~ 1 + logCS,
  ~ 1 + logCS
)

CAC_2breaks_int_mcp <- mcp(CAC_2breaks_int_model, LM_MCP_array)
plot(CAC_2breaks_int_mcp)
summary(CAC_2breaks_int_mcp)
#

#Null model - no break points and single slope/intercept
Null_CAC_eco_model = list(
  CACs ~ (1|simple_ecomorph) + logCS
)

Null_CAC_eco_mcp <- mcp(Null_CAC_eco_model, LM_MCP_array)
plot(Null_CAC_eco_mcp)
summary(Null_CAC_LM_mcp)
ranef(Null_CAC_eco_mcp)
#

Null_CAC_LM_mcp$loo <- loo(Null_CAC_LM_mcp)
Null_CAC_eco_mcp$loo <- loo(Null_CAC_eco_mcp)
CAC_LM_mcp$loo <- loo(CAC_LM_mcp)
CAC_LM_int_mcp$loo <- loo(CAC_LM_int_mcp)
CAC_LM_sp_mcp$loo <- loo(CAC_LM_sp_mcp)
CAC_LM_sp_mcp_V2$loo <- loo(CAC_LM_sp_mcp_V2)
CAC_LM_eco_mcp_V3$loo <- loo(CAC_LM_eco_mcp_V3)
CAC_LM_eco_mcp$loo <- loo(CAC_LM_eco_mcp)
CAC_LM_eco_mcp_V2$loo <- loo(CAC_LM_eco_mcp_V2)
CAC_2breaks_int_mcp$loo <- loo(CAC_2breaks_int_mcp)


loo::loo_compare(Null_CAC_LM_mcp$loo, Null_CAC_eco_mcp$loo, CAC_LM_mcp$loo, CAC_LM_int_mcp$loo, CAC_LM_sp_mcp$loo, CAC_LM_sp_mcp_V2$loo, CAC_LM_eco_mcp$loo, CAC_LM_eco_mcp_V2$loo)

loo::loo_compare(Null_CAC_LM_mcp$loo, CAC_LM_mcp$loo, CAC_LM_int_mcp$loo, CAC_LM_sp_mcp$loo, CAC_LM_eco_mcp$loo, CAC_LM_log_eco_mcp$loo)
##Simple and complex ecomorph differences in slope and intercept are much better fits than all other models. This is not the same as the skull table curve alone!

loo::loo_compare(Null_CAC_LM_mcp$loo, CAC_LM_int_mcp$loo, CAC_LM_eco_mcp$loo, CAC_LM_eco_mcp_V3$loo)

loo::loo_compare(Null_CAC_LM_mcp$loo, CAC_LM_int_mcp$loo, CAC_LM_sp_mcp$loo, CAC_LM_eco_mcp$loo)
###

###Dorsal Face
PCA <- DorsalFace_20_PCA
CACs <- DorsalFace_20_CACs

LM_MCP_array <- data.frame(PCA$x,CACs, log(Croc_Lateral_GPA$Csize), lateral_covariates$species, lateral_covariates$Ecomorph, lateral_covariates$simple_ecomorph, stringsAsFactors=FALSE)
dimnames(LM_MCP_array)[[2]][42] <- "logCS"
dimnames(LM_MCP_array)[[2]][43] <- "species"
dimnames(LM_MCP_array)[[2]][44] <- "ecomorph"
dimnames(LM_MCP_array)[[2]][45] <- "simple_ecomorph"

#Null model - no break points and single slope/intercept
Null_CAC_model = list(
  CACs ~ logCS
)

Null_CAC_LM_mcp <- mcp(Null_CAC_model, LM_MCP_array)
plot(Null_CAC_LM_mcp)
#

#model for single break point with multiple slopes
CAC_LM_model = list(
  CACs ~ 1 + logCS,
  ~ 0 + logCS
  
)

poly(LM_MCP_array[c("CACs","logCS")],degree=3)

CAC_LM_mcp <- mcp(CAC_LM_model, LM_MCP_array)
plot(CAC_LM_mcp)
summary(CAC_LM_mcp)
#

#model for single break point with multiple intercepts and slopes
CAC_LM_int_model = list(
  CACs ~ 1 + logCS,
  ~ 1 + logCS
)

CAC_LM_int_mcp <- mcp(CAC_LM_int_model, LM_MCP_array)
plot(CAC_LM_int_mcp)
summary(CAC_LM_int_mcp)
#

#model for single break point with variance among species and multiple intercepts and slopes in both lines
CAC_LM_sp_model_V2 = list(
  CACs ~ 1 + (1|species) + logCS,
  ~ (0|species) + logCS
)

CAC_LM_sp_mcp_V2 <- mcp(CAC_LM_sp_model_V2, LM_MCP_array)
plot(CAC_LM_sp_mcp_V2)
summary(CAC_LM_sp_mcp_V2)
#

#model for single break point with variance among species with multiple intercepts and slopes for only the second line
CAC_LM_sp_model = list(
  CACs ~ 1 + logCS,
  1 + (1|species) ~ (0|species) + logCS
)

CAC_LM_sp_mcp <- mcp(CAC_LM_sp_model, LM_MCP_array)
plot(CAC_LM_sp_mcp)
plot(CAC_LM_sp_mcp, facet_by = "species")
summary(CAC_LM_sp_mcp)
#

#model for single break point with variance among species and multiple intercepts and slopes in both lines
CAC_LM_eco_model_V2 = list(
  CACs ~ 1 + (1|simple_ecomorph) + logCS,
  ~ (0|simple_ecomorph) + logCS
)

CAC_LM_eco_mcp_V2 <- mcp(CAC_LM_eco_model_V2, LM_MCP_array)
plot(CAC_LM_eco_mcp_V2)
summary(CAC_LM_eco_mcp_V2)

#model for single break point with variance among species and multiple intercepts and slopes in the second line only
CAC_LM_eco_model = list(
  CACs ~ 1 + logCS,
  1 + (1|simple_ecomorph) ~ (0|simple_ecomorph) + logCS
)

CAC_LM_eco_mcp <- mcp(CAC_LM_eco_model, LM_MCP_array)
plot(CAC_LM_eco_mcp)
plot(CAC_LM_eco_mcp, facet_by = "simple_ecomorph")
summary(CAC_LM_eco_mcp)
ranef(CAC_LM_eco_mcp)

#model for two break points with multiple intercepts and slopes
CAC_2breaks_int_model = list(
  CACs ~ 2 + logCS,
  ~ 1 + logCS,
  ~ 1 + logCS
)

CAC_2breaks_int_mcp <- mcp(CAC_2breaks_int_model, LM_MCP_array)
plot(CAC_2breaks_int_mcp)
ummary(CAC_2breaks_int_mcp)
#

#Null model - no break points and single slope/intercept
Null_CAC_eco_model = list(
  CACs ~ (1|simple_ecomorph) + logCS
)

Null_CAC_eco_mcp <- mcp(Null_CAC_eco_model, LM_MCP_array)
plot(Null_CAC_eco_mcp)
#

Null_CAC_LM_mcp$loo <- loo(Null_CAC_LM_mcp)
Null_CAC_eco_mcp$loo <- loo(Null_CAC_eco_mcp)
CAC_LM_mcp$loo <- loo(CAC_LM_mcp)
CAC_LM_int_mcp$loo <- loo(CAC_LM_int_mcp)
CAC_LM_sp_mcp$loo <- loo(CAC_LM_sp_mcp)
CAC_LM_sp_mcp_V2$loo <- loo(CAC_LM_sp_mcp_V2)
CAC_LM_eco_mcp$loo <- loo(CAC_LM_eco_mcp)
CAC_LM_eco_mcp_V2$loo <- loo(CAC_LM_eco_mcp_V2)


loo::loo_compare(Null_CAC_LM_mcp$loo, CAC_LM_mcp$loo, CAC_LM_int_mcp$loo, CAC_LM_sp_mcp$loo, CAC_LM_sp_mcp_V2$loo, CAC_LM_eco_mcp$loo, CAC_LM_eco_mcp_V2$loo, Null_CAC_eco_mcp$loo)

loo::loo_compare(Null_CAC_LM_mcp$loo, CAC_LM_mcp$loo, CAC_LM_int_mcp$loo, CAC_LM_sp_mcp$loo, CAC_LM_eco_mcp$loo)
###

### Skull Table Curve alone
PCA <- SkullTable_20_PCA
CACs <- SkullTable_20_CACs

MCP_array <- data.frame(PCA$x, CACs, log(Croc_Lateral_GPA$Csize), lateral_covariates$species,lateral_covariates$Ecomorph, lateral_covariates$simple_ecomorph, stringsAsFactors=FALSE)
dimnames(MCP_array)[[2]][42] <- "logCS"
dimnames(MCP_array)[[2]][43] <- "species"
dimnames(MCP_array)[[2]][44] <- "ecomorph"
dimnames(MCP_array)[[2]][45] <- "simple_ecomorph"

#Null model - no break points and single slope/intercept
Null_CAC_model = list(
  CACs ~ logCS
)

Null_CAC_mcp <- mcp(Null_CAC_model, MCP_array)
plot(Null_CAC_mcp)
#

#model for single break point with multiple slopes
CACs_model = list(
  CACs ~ 1 + logCS,
  ~ 0 + logCS
)

CACs_mcp <- mcp(CACs_model, MCP_array)
plot(CACs_mcp)
summary(CACs_mcp)
#

#model for single break point with multiple intercepts and slopes
CACs_int_model = list(
  CACs ~ 1 + logCS,
  ~ 1 + logCS
)

CACs_int_mcp <- mcp(CACs_int_model, MCP_array)
plot(CACs_int_mcp)
summary(CACs_int_mcp)
#

#model for single break point with variance among species and multiple intercepts and slopes
CACs_sp_model = list(
  CACs ~ 1 + logCS,
  1 + (1|species) ~ (0|species) + logCS
)

CACs_sp_mcp <- mcp(CACs_sp_model, MCP_array)
plot(CACs_sp_mcp)
plot(CACs_sp_mcp, facet_by = "species")
summary(CACs_sp_mcp)
#

#model for single break point with variance among species and multiple intercepts and slopes
CACs_eco_model = list(
  CACs ~ 1 + logCS,
  1 + (1|simple_ecomorph) ~ (0|simple_ecomorph) + logCS
)

CACs_eco_mcp <- mcp(CACs_eco_model, MCP_array)
plot(CACs_eco_mcp)
plot(CACs_eco_mcp, facet_by = "simple_ecomorph")
summary(CACs_eco_mcp)
ranef(CACs_eco_mcp)
#


Null_CAC_mcp$loo <- loo(Null_CAC_mcp)
CACs_mcp$loo <- loo(CACs_mcp)
CACs_int_mcp$loo <- loo(CACs_int_mcp)
CACs_sp_mcp$loo <- loo(CACs_sp_mcp)
CACs_eco_mcp$loo <- loo(CACs_eco_mcp)


loo::loo_compare(Null_CAC_mcp$loo, CACs_int_mcp$loo, CACs_sp_mcp$loo, CACs_eco_mcp$loo)
###