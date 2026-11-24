library(Hmsc)

da = readRDS("outputs/sim_survey_100S_5V.rds")
meta = da$meta
studyDesign = data.frame(site = as.factor(da$visits$cell_id), doy = as.factor(da$visits$doy))
sites = unique(da$visits$cell_id)
xy = as.matrix(da$visits[match(sites,da$visits$cell_id),c("col","row")])
colnames(xy) = c("x-coordinate","y-coordinate")
rownames(xy) = sites
XData = da$cells[match(da$visits$cell_id,da$cells$cell_id),c("forest_prop_z","edge_density_z","forest_buf1km_z")]
XData$doy = da$visits$doy
ny = nrow(XData)

spp = unique(da$detections$species_id)
ns = length(spp)
Y = matrix(0,nrow = ny, ncol = ns)
colnames(Y) = spp
for(i in 1:ny){
  sel = which((studyDesign$site[i]==da$detections$cell_id) &
                (studyDesign$doy[i]==da$detections$doy))
  da1 = da$detections[sel,]
  Y[i,da1$species_id[da1$max_score>=0.9]] = 1
}

XFormula = ~forest_prop_z+edge_density_z+forest_buf1km_z+poly(doy,degree = 2,raw = T)
rL.site = HmscRandomLevel(sData = xy)

m = Hmsc(XData = XData, XFormula = XFormula,
         Y = Y, studyDesign = studyDesign,
         ranLevels = list(site = rL.site),
         distr = "probit")
models = list(pa = m)
save(models,file = "models/unfitted_models.RData")
