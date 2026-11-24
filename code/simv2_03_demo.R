# simv2_08_demo.R
# Example session. Run from the folder that holds the simv2_0*.R files.

for (f in sort(list.files(".", pattern = "^simv2_0[0-2].*\\.R$"))) source(f)

cfg <- sim_config(seed = 1)
#  Overrides to test diff settings of cfg
#   sim_config(occupancy = list(beta_mean = c(forest = 1, edge = 0, buffer = 0.5)))
#   sim_config(community = list(n_visitor_only = 5))     # bring visitor-only species back (extra to the 50)
#   sim_config(occupancy = list(field_scale_m = c(NA, 300, 1000, 3000, 6000)))   # 6 km instead of 10 km
#   sim_config(occupancy = list(field_balanced = TRUE))  # near-equal number of species per scale
#   sim_config(occupancy = list(p_unaffected = 0))       # switch off the fragmentation draw: every species gets betas
#   sim_config(occupancy = list(p_unaffected = 0.35))    # more species unaffected by fragmentation
#   sim_config(covariate_file = "path/to/other_grid.csv")

world <- sim_init(cfg)

sp <- sim_truth()$species
sim_write_species_table(sp, "outputs/simulated_species.csv")
#getwd()    # shows the folder where species.csv was saved

# One 24 h recording at row 120, column 130 on day 140. Only max_score is returned.
simulate(120, 130, 140)

# A small design: three cells visited on three days.
# Water cells cannot be sampled, so pick cells from the non-water ones.
land_cells <- which(!world$land$water, arr.ind = TRUE)   # columns: row (= i), col (= j)
pick <- land_cells[c(1000, 20000, 40000), ]
design <- data.frame(i = pick[, "row"], j = pick[, "col"], doy = c(140, 150, 180))
head(sim_collect(design))

# Truth (Don't given to the sampler)
sp <- sim_truth()$species
head(sp)
sim_plot_truth(sp$species_id[1])
sim_truth_day(120, 130, 140)[1:10, ]

# Edit the community if needed
# sim_write_species_table(sp, "species.csv")   # edit in a spreadsheet, then:
# world2 <- sim_init(cfg, species = sim_read_species_table("species.csv"))

# Checks
v <- sim_validate()
