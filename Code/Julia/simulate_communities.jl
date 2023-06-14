################################################
## Prepare stress vector ##
################################################
ctrl_harvest = 3.0 #follows Barauh et al 2021 (typical collapse at 1.5-2.8)
y_harvest = [fill(0.0,201);(0:(ctrl_harvest-0)/99:ctrl_harvest)] #control increases from t200 to t300
signal_harvest = LinearInterpolation((0:1:300), y_harvest)

ctrl_invasive = 5.0 #follows Dakos 2018 (typical collapse at 2.5-3.0)
y_invasive = [fill(0.0,201);(0:(ctrl_invasive-0)/99:ctrl_invasive)] #control increases from t200 to t300
signal_invasive = LinearInterpolation((0:1:300), y_invasive)

y_null = fill(0.0,301) #control increases from t200 to t300
signal_null = LinearInterpolation((0:1:300), y_null)

motifs = (1,2,6,9,10)

################################################
## 5 Species ##
################################################
#Harvesting model auto
out_df =DataFrame([[],[],[],[],[]], [:community,:collapse,:stress_param,:max_eigval,:motif]) #array to store jacobian information

for i in motifs
    
    #i = 1
    motif_dat = RData.load(string("Data/networks/5_spp/web_",i,".RData"),convert=true)["out_file"] #load motif networks generated in R
    #motif_jac_dat = extract_max_eigvec(motif_dat,0.0:0.01:3.0, 100, "auto"; model = "harvest",abs_max = false,concrete_jac = true) #simulate deterministic model for each community
    motif_jac_dat = extract_max_eigvec_alt(motif_dat[1:100],0.0:0.01:ctrl_harvest, 100, "auto"; model = "harvest",abs_max = false) #simulate deterministic model for each community
    motif_jac_dat.motif .= string(i) #label motif
    motif_sim_dat = motif_sim(motif_dat[1:100],missing,"auto",repeat([0.1],length(motif_dat[i]["Neq"])),(0,300),signal_harvest,100,noise_col = "white",saveat=1.0)  #simulate stochastic communities
    #in order, arguments are: parameter data (interaction matrix, starting abundances, carrying capacities etc), starting abundance (if "missing", uses data from the previous argument),
    #species to stress (if "auto", targets the most abundant species in the third trophic level), noise scalar, start-endpoint of timeseries,
    #stress vector, number of simulations (per community), noise colour (white or correlated), saveat what resolution (dt = 0.1) 
    # motif_sim_dat = motif_sim(motif_dat,harvested_spp_ls[string(i)],[0.1,0.1,0.1,0.1],(0,500),signal_motif,100,noise_col = "white",saveat=1.0) 
    motif_csv = extract_ts_for_csv(motif_sim_dat,string(i),100:300,5) #crop to t100 - t300, label as motif "i" and round to 5 decimal places
    CSV.write(string("Data/simulations/5_spp/harvest/stress/motif_",i,".csv"),motif_csv,compress=true) #save simulations
    append!(out_df,motif_jac_dat) #add jacobian info to expanding dataframe

    motif_sim_unstress_dat = motif_sim(motif_dat[1:100],missing,"auto",repeat([0.1],length(motif_dat[i]["Neq"])),(0,300),signal_null,100,noise_col = "white",saveat=1.0)  #simulate stochastic communities but unstressed
    motif_unstress_csv = extract_ts_for_csv(motif_sim_unstress_dat,string(i),100:300,5) 
    CSV.write(string("Data/simulations/5_spp/harvest/unstressed/motif_",i,".csv"),motif_unstress_csv,compress=false) #save simulations

if i == last(motifs) 
    CSV.write(string("Data/simulations/5_spp/harvest/stress/jacobian_5_spp.csv"),out_df,compress=false) #save out jacobian info on last iteration
end
end

#Harvesting model 20%
out_df =DataFrame([[],[],[],[],[]], [:community,:collapse,:stress_param,:max_eigval,:motif]) #array to store jacobian information
harv_spp = [4]
for i in motifs
    
    #i = 1
    motif_dat = RData.load(string("Data/networks/5_spp/web_",i,".RData"),convert=true)["out_file"] #load motif networks generated in R
    #motif_jac_dat = extract_max_eigvec(motif_dat,0.0:0.01:3.0, 100, "auto"; model = "harvest",abs_max = false,concrete_jac = true) #simulate deterministic model for each community
    motif_jac_dat = extract_max_eigvec_alt(motif_dat[1:100],0.0:0.01:ctrl_harvest, 100, harv_spp; model = "harvest",abs_max = false) #simulate deterministic model for each community
    motif_jac_dat.motif .= string(i) #label motif
    motif_sim_dat = motif_sim(motif_dat[1:100],missing,harv_spp,repeat([0.1],length(motif_dat[i]["Neq"])),(0,300),signal_harvest,100,noise_col = "white",saveat=1.0)  #simulate stochastic communities
    #in order, arguments are: parameter data (interaction matrix, starting abundances, carrying capacities etc), starting abundance (if "missing", uses data from the previous argument),
    #species to stress (if "auto", targets the most abundant species in the third trophic level), noise scalar, start-endpoint of timeseries,
    #stress vector, number of simulations (per community), noise colour (white or correlated), saveat what resolution (dt = 0.1) 
    # motif_sim_dat = motif_sim(motif_dat,harvested_spp_ls[string(i)],[0.1,0.1,0.1,0.1],(0,500),signal_motif,100,noise_col = "white",saveat=1.0) 
    motif_csv = extract_ts_for_csv(motif_sim_dat,string(i),100:300,5) #crop to t100 - t300, label as motif "i" and round to 5 decimal places
    CSV.write(string("Data/simulations/5_spp/harvest20/stress/motif_",i,"_20percent.csv"),motif_csv,compress=true) #save simulations
    append!(out_df,motif_jac_dat) #add jacobian info to expanding dataframe

    motif_sim_unstress_dat = motif_sim(motif_dat[1:100],missing,"auto",repeat([0.1],length(motif_dat[i]["Neq"])),(0,300),signal_null,100,noise_col = "white",saveat=1.0)  #simulate stochastic communities but unstressed
    motif_unstress_csv = extract_ts_for_csv(motif_sim_unstress_dat,string(i),100:300,5) 
    CSV.write(string("Data/simulations/5_spp/harvest20/unstressed/motif_",i,"_20percent.csv"),motif_unstress_csv,compress=false) #save simulations

if i == last(motifs) 
    CSV.write(string("Data/simulations/5_spp/harvest20/stress/jacobian_5_spp_20percent.csv"),out_df,compress=false) #save out jacobian info on last iteration
end
end

#Invasive
out_df =DataFrame([[],[],[],[],[]], [:community,:collapse,:stress_param,:max_eigval,:motif]) #array to store jacobian information

for i in motifs
    motif_dat = RData.load(string("Data/networks/5_spp/web_",i,".RData"),convert=true)["out_file"] #load motif networks generated in R

    inv_motif_dat = Vector{DictoVec{}}(undef,length(motif_dat[1:100])) #invasive model requires extra community info. Prepare array to hold this info
for j in eachindex(motif_dat[1:100])

    inv_motif_dat[j] = RData.DictoVec([-1*motif_dat[j]["A_matrix"],
    rand(Distributions.Uniform(0.9,1.1),length(motif_dat[j]["Neq"])),
    rand(Distributions.Uniform(5,15),length(motif_dat[j]["Neq"])),
    repeat([0.1],length(motif_dat[j]["Neq"])),
    #rand(Distributions.Uniform(0,1),length(motif_dat[j]["Neq"])), #random eta
    repeat([1.0],length(motif_dat[j]["Neq"])), #set eta
    motif_dat[j]["Neq"]*10,
    motif_dat[j]["tlvl"]],
    ["A_matrix","R","K","imm","eta","Neq","tlvl"]) #generate additional random info. A_matrix = loaded interaction matrix, R = growth rate,
    #K = carrying capacity, imm = immigration term, eta = stress scalar, Neq = starting abundance, tlvl = trophic level of each species
end

    motif_jac_dat = extract_max_eigvec_alt(inv_motif_dat,0.0:0.01:ctrl_invasive, 100, "auto"; model = "invasive",abs_max = false) #simulate deterministic model for each community
    motif_jac_dat.motif .= string(i) #label motif
    motif_sim_dat = invasive_sim(inv_motif_dat,missing,"auto",repeat([0.1],length(inv_motif_dat[i]["Neq"])),(0,300),signal_invasive,100,noise_col = "white",saveat=1.0) #simulate stochastic communities.
    #only difference to above is that noise scalar is set generically to be length of number of species
    motif_csv = extract_ts_for_csv(motif_sim_dat,string(i),100:300,5) #crop to t100 - t400 and label as motif "1"
    CSV.write(string("Data/simulations/5_spp/invasive/stress/motif_",i,".csv"),motif_csv,compress=true) #save simulations
    append!(out_df,motif_jac_dat)

    motif_sim_unstress_dat = invasive_sim(inv_motif_dat,missing,"auto",repeat([0.1],length(inv_motif_dat[i]["Neq"])),(0,300),signal_null,100,noise_col = "white",saveat=1.0) #simulate stochastic communities.
    motif_unstress_csv = extract_ts_for_csv(motif_sim_unstress_dat,string(i),100:300,5) #crop to t100 - t400 and label as motif "1"
    CSV.write(string("Data/simulations/5_spp/invasive/unstressed/motif_",i,".csv"),motif_unstress_csv,compress=true) #save simulations

if i == last(motifs)  
    CSV.write(string("Data/simulations/5_spp/invasive/stress/jacobian_5_spp.csv"),out_df,compress=false) #save out jacobian info on last iteration
end
end

################################################
## 10 Species ##
################################################
#Harvesting model auto
out_df =DataFrame([[],[],[],[],[]], [:community,:collapse,:stress_param,:max_eigval,:motif]) #array to store jacobian information

for i in motifs
    #i = 1
    motif_dat = RData.load(string("Data/networks/10_spp/web_",i,".RData"),convert=true)["out_file"] #load motif networks generated in R
    #motif_jac_dat = extract_max_eigvec(motif_dat,0.0:0.01:3.0, 100, "auto"; model = "harvest",abs_max = false,concrete_jac = true) #simulate deterministic model for each community
    motif_jac_dat = extract_max_eigvec_alt(motif_dat[1:100],0.0:0.01:ctrl_harvest, 100, "auto"; model = "harvest",abs_max = false) #simulate deterministic model for each community
    motif_jac_dat.motif .= string(i) #label motif
    motif_sim_dat = motif_sim(motif_dat[1:100],missing,"auto",repeat([0.1],length(motif_dat[i]["Neq"])),(0,300),signal_harvest,100,noise_col = "white",saveat=1.0)  #simulate stochastic communities
    #in order, arguments are: parameter data (interaction matrix, starting abundances, carrying capacities etc), starting abundance (if "missing", uses data from the previous argument),
    #species to stress (if "auto", targets the most abundant species in the third trophic level), noise scalar, start-endpoint of timeseries,
    #stress vector, number of simulations (per community), noise colour (white or correlated), saveat what resolution (dt = 0.1) 
    # motif_sim_dat = motif_sim(motif_dat,harvested_spp_ls[string(i)],[0.1,0.1,0.1,0.1],(0,500),signal_motif,100,noise_col = "white",saveat=1.0) 
    motif_csv = extract_ts_for_csv(motif_sim_dat,string(i),100:300,5) #crop to t100 - t300, label as motif "i" and round to 5 decimal places
    CSV.write(string("Data/simulations/10_spp/harvest/stress/motif_",i,".csv"),motif_csv,compress=true) #save simulations
    append!(out_df,motif_jac_dat) #add jacobian info to expanding dataframe

    motif_sim_unstress_dat = motif_sim(motif_dat[1:100],missing,"auto",repeat([0.1],length(motif_dat[i]["Neq"])),(0,300),signal_null,100,noise_col = "white",saveat=1.0)  #simulate stochastic communities but unstressed
    motif_unstress_csv = extract_ts_for_csv(motif_sim_unstress_dat,string(i),100:300,5) 
    CSV.write(string("Data/simulations/10_spp/harvest/unstressed/motif_",i,".csv"),motif_unstress_csv,compress=true) #save simulations

if i == last(motifs)   
    CSV.write(string("Data/simulations/10_spp/harvest/stress/jacobian_10_spp.csv"),out_df,compress=false) #save out jacobian info on last iteration
end
end

#Harvesting model 20%
out_df =DataFrame([[],[],[],[],[]], [:community,:collapse,:stress_param,:max_eigval,:motif]) #array to store jacobian information
harv_spp = [8,9]

for i in motifs
    #i = 1
    motif_dat = RData.load(string("Data/networks/10_spp/web_",i,".RData"),convert=true)["out_file"] #load motif networks generated in R
    #motif_jac_dat = extract_max_eigvec(motif_dat,0.0:0.01:3.0, 100, "auto"; model = "harvest",abs_max = false,concrete_jac = true) #simulate deterministic model for each community
    motif_jac_dat = extract_max_eigvec_alt(motif_dat[1:100],0.0:0.01:ctrl_harvest, 100, harv_spp; model = "harvest",abs_max = false) #simulate deterministic model for each community
    motif_jac_dat.motif .= string(i) #label motif
    motif_sim_dat = motif_sim(motif_dat[1:100],missing,harv_spp,repeat([0.1],length(motif_dat[i]["Neq"])),(0,300),signal_harvest,100,noise_col = "white",saveat=1.0)  #simulate stochastic communities
    #in order, arguments are: parameter data (interaction matrix, starting abundances, carrying capacities etc), starting abundance (if "missing", uses data from the previous argument),
    #species to stress (if "auto", targets the most abundant species in the third trophic level), noise scalar, start-endpoint of timeseries,
    #stress vector, number of simulations (per community), noise colour (white or correlated), saveat what resolution (dt = 0.1) 
    # motif_sim_dat = motif_sim(motif_dat,harvested_spp_ls[string(i)],[0.1,0.1,0.1,0.1],(0,500),signal_motif,100,noise_col = "white",saveat=1.0) 
    motif_csv = extract_ts_for_csv(motif_sim_dat,string(i),100:300,5) #crop to t100 - t300, label as motif "i" and round to 5 decimal places
    CSV.write(string("Data/simulations/10_spp/harvest20/stress/motif_",i,"_20percent.csv"),motif_csv,compress=true) #save simulations
    append!(out_df,motif_jac_dat) #add jacobian info to expanding dataframe

    motif_sim_unstress_dat = motif_sim(motif_dat[1:100],missing,"auto",repeat([0.1],length(motif_dat[i]["Neq"])),(0,300),signal_null,100,noise_col = "white",saveat=1.0)  #simulate stochastic communities but unstressed
    motif_unstress_csv = extract_ts_for_csv(motif_sim_unstress_dat,string(i),100:300,5) 
    CSV.write(string("Data/simulations/10_spp/harvest20/unstressed/motif_",i,"_20percent.csv"),motif_unstress_csv,compress=true) #save simulations

if i == last(motifs)   
    CSV.write(string("Data/simulations/10_spp/harvest20/stress/jacobian_10_spp_20percent.csv"),out_df,compress=false) #save out jacobian info on last iteration
end
end

#Invasive
out_df =DataFrame([[],[],[],[],[]], [:community,:collapse,:stress_param,:max_eigval,:motif]) #array to store jacobian information

for i in motifs
    motif_dat = RData.load(string("Data/networks/10_spp/web_",i,".RData"),convert=true)["out_file"] #load motif networks generated in R

    inv_motif_dat = Vector{DictoVec{}}(undef,length(motif_dat[1:100])) #invasive model requires extra community info. Prepare array to hold this info
for j in eachindex(motif_dat[1:100])

    inv_motif_dat[j] = RData.DictoVec([-1*motif_dat[j]["A_matrix"],
    rand(Distributions.Uniform(0.9,1.1),length(motif_dat[j]["Neq"])),
    rand(Distributions.Uniform(5,15),length(motif_dat[j]["Neq"])),
    repeat([0.1],length(motif_dat[j]["Neq"])),
    #rand(Distributions.Uniform(0,1),length(motif_dat[j]["Neq"])), #random eta
    repeat([1.0],length(motif_dat[j]["Neq"])), #set eta
    motif_dat[j]["Neq"]*10,
    motif_dat[j]["tlvl"]],
    ["A_matrix","R","K","imm","eta","Neq","tlvl"]) #generate additional random info. A_matrix = loaded interaction matrix, R = growth rate,
    #K = carrying capacity, imm = immigration term, eta = stress scalar, Neq = starting abundance, tlvl = trophic level of each species
end

    motif_jac_dat = extract_max_eigvec_alt(inv_motif_dat,0.0:0.01:ctrl_invasive, 100, "auto"; model = "invasive",abs_max = false) #simulate deterministic model for each community
    motif_jac_dat.motif .= string(i) #label motif
    motif_sim_dat = invasive_sim(inv_motif_dat,missing,"auto",repeat([0.1],length(inv_motif_dat[i]["Neq"])),(0,300),signal_invasive,100,noise_col = "white",saveat=1.0) #simulate stochastic communities.
    #only difference to above is that noise scalar is set generically to be length of number of species
    motif_csv = extract_ts_for_csv(motif_sim_dat,string(i),100:300,5) #crop to t100 - t400 and label as motif "1"
    CSV.write(string("Data/simulations/10_spp/invasive/stress/motif_",i,".csv"),motif_csv,compress=true) #save simulations
    append!(out_df,motif_jac_dat)

    motif_sim_unstress_dat = invasive_sim(inv_motif_dat,missing,"auto",repeat([0.1],length(inv_motif_dat[i]["Neq"])),(0,300),signal_null,100,noise_col = "white",saveat=1.0) #simulate stochastic communities.
    motif_unstress_csv = extract_ts_for_csv(motif_sim_unstress_dat,string(i),100:300,5) #crop to t100 - t400 and label as motif "1"
    CSV.write(string("Data/simulations/10_spp/invasive/unstressed/motif_",i,".csv"),motif_unstress_csv,compress=true) #save simulations

if i == last(motifs)   
    CSV.write(string("Data/simulations/10_spp/invasive/stress/jacobian_10_spp.csv"),out_df,compress=false) #save out jacobian info on last iteration
end
end

################################################
## 15 Species ##
################################################
#Harvesting model auto
out_df =DataFrame([[],[],[],[],[]], [:community,:collapse,:stress_param,:max_eigval,:motif]) #array to store jacobian information

for i in motifs
    #i = 1
    motif_dat = RData.load(string("Data/networks/15_spp/web_",i,".RData"),convert=true)["out_file"] #load motif networks generated in R
    #motif_jac_dat = extract_max_eigvec(motif_dat,0.0:0.01:3.0, 100, "auto"; model = "harvest",abs_max = false,concrete_jac = true) #simulate deterministic model for each community
    motif_jac_dat = extract_max_eigvec_alt(motif_dat[1:100],0.0:0.01:ctrl_harvest, 100, "auto"; model = "harvest",abs_max = false) #simulate deterministic model for each community
    motif_jac_dat.motif .= string(i) #label motif
    motif_sim_dat = motif_sim(motif_dat[1:100],missing,"auto",repeat([0.1],length(motif_dat[i]["Neq"])),(0,300),signal_harvest,100,noise_col = "white",saveat=1.0)  #simulate stochastic communities
    #in order, arguments are: parameter data (interaction matrix, starting abundances, carrying capacities etc), starting abundance (if "missing", uses data from the previous argument),
    #species to stress (if "auto", targets the most abundant species in the third trophic level), noise scalar, start-endpoint of timeseries,
    #stress vector, number of simulations (per community), noise colour (white or correlated), saveat what resolution (dt = 0.1) 
    # motif_sim_dat = motif_sim(motif_dat,harvested_spp_ls[string(i)],[0.1,0.1,0.1,0.1],(0,500),signal_motif,100,noise_col = "white",saveat=1.0) 
    motif_csv = extract_ts_for_csv(motif_sim_dat,string(i),100:300,5) #crop to t100 - t300, label as motif "i" and round to 5 decimal places
    CSV.write(string("Data/simulations/15_spp/harvest/stress/motif_",i,".csv"),motif_csv,compress=true) #save simulations
    append!(out_df,motif_jac_dat) #add jacobian info to expanding dataframe

    motif_sim_unstress_dat = motif_sim(motif_dat[1:100],missing,"auto",repeat([0.1],length(motif_dat[i]["Neq"])),(0,300),signal_null,100,noise_col = "white",saveat=1.0)  #simulate stochastic communities but unstressed
    motif_unstress_csv = extract_ts_for_csv(motif_sim_unstress_dat,string(i),100:300,5) 
    CSV.write(string("Data/simulations/15_spp/harvest/unstressed/motif_",i,".csv"),motif_unstress_csv,compress=true) #save simulations

if i == last(motifs)   
    CSV.write(string("Data/simulations/15_spp/harvest/stress/jacobian_15_spp.csv"),out_df,compress=false) #save out jacobian info on last iteration
end
end

#Harvesting model 20%
out_df =DataFrame([[],[],[],[],[]], [:community,:collapse,:stress_param,:max_eigval,:motif]) #array to store jacobian information
harv_spp = [11,12,13]

for i in motifs
    #i = 1
    motif_dat = RData.load(string("Data/networks/15_spp/web_",i,".RData"),convert=true)["out_file"] #load motif networks generated in R
    #motif_jac_dat = extract_max_eigvec(motif_dat,0.0:0.01:3.0, 100, "auto"; model = "harvest",abs_max = false,concrete_jac = true) #simulate deterministic model for each community
    motif_jac_dat = extract_max_eigvec_alt(motif_dat[1:100],0.0:0.01:ctrl_harvest, 100, harv_spp; model = "harvest",abs_max = false) #simulate deterministic model for each community
    motif_jac_dat.motif .= string(i) #label motif
    motif_sim_dat = motif_sim(motif_dat[1:100],missing,harv_spp,repeat([0.1],length(motif_dat[i]["Neq"])),(0,300),signal_harvest,100,noise_col = "white",saveat=1.0)  #simulate stochastic communities
    #in order, arguments are: parameter data (interaction matrix, starting abundances, carrying capacities etc), starting abundance (if "missing", uses data from the previous argument),
    #species to stress (if "auto", targets the most abundant species in the third trophic level), noise scalar, start-endpoint of timeseries,
    #stress vector, number of simulations (per community), noise colour (white or correlated), saveat what resolution (dt = 0.1) 
    # motif_sim_dat = motif_sim(motif_dat,harvested_spp_ls[string(i)],[0.1,0.1,0.1,0.1],(0,500),signal_motif,100,noise_col = "white",saveat=1.0) 
    motif_csv = extract_ts_for_csv(motif_sim_dat,string(i),100:300,5) #crop to t100 - t300, label as motif "i" and round to 5 decimal places
    CSV.write(string("Data/simulations/15_spp/harvest20/stress/motif_",i,"_20percent.csv"),motif_csv,compress=true) #save simulations
    append!(out_df,motif_jac_dat) #add jacobian info to expanding dataframe

    motif_sim_unstress_dat = motif_sim(motif_dat[1:100],missing,"auto",repeat([0.1],length(motif_dat[i]["Neq"])),(0,300),signal_null,100,noise_col = "white",saveat=1.0)  #simulate stochastic communities but unstressed
    motif_unstress_csv = extract_ts_for_csv(motif_sim_unstress_dat,string(i),100:300,5) 
    CSV.write(string("Data/simulations/15_spp/harvest20/unstressed/motif_",i,"_20percent.csv"),motif_unstress_csv,compress=true) #save simulations

if i == last(motifs)   
    CSV.write(string("Data/simulations/15_spp/harvest20/stress/jacobian_15_spp_20percent.csv"),out_df,compress=false) #save out jacobian info on last iteration
end
end

#Invasive
out_df =DataFrame([[],[],[],[],[]], [:community,:collapse,:stress_param,:max_eigval,:motif]) #array to store jacobian information

for i in motifs
    motif_dat = RData.load(string("Data/networks/15_spp/web_",i,".RData"),convert=true)["out_file"] #load motif networks generated in R

    inv_motif_dat = Vector{DictoVec{}}(undef,length(motif_dat[1:100])) #invasive model requires extra community info. Prepare array to hold this info
for j in eachindex(motif_dat[1:100])

    inv_motif_dat[j] = RData.DictoVec([-1*motif_dat[j]["A_matrix"],
    rand(Distributions.Uniform(0.9,1.1),length(motif_dat[j]["Neq"])),
    rand(Distributions.Uniform(5,15),length(motif_dat[j]["Neq"])),
    repeat([0.1],length(motif_dat[j]["Neq"])),
    #rand(Distributions.Uniform(0,1),length(motif_dat[j]["Neq"])), #random eta
    repeat([1.0],length(motif_dat[j]["Neq"])), #set eta
    motif_dat[j]["Neq"]*10,
    motif_dat[j]["tlvl"]],
    ["A_matrix","R","K","imm","eta","Neq","tlvl"]) #generate additional random info. A_matrix = loaded interaction matrix, R = growth rate,
    #K = carrying capacity, imm = immigration term, eta = stress scalar, Neq = starting abundance, tlvl = trophic level of each species
end

    motif_jac_dat = extract_max_eigvec_alt(inv_motif_dat,0.0:0.01:ctrl_invasive, 100, "auto"; model = "invasive",abs_max = false) #simulate deterministic model for each community
    motif_jac_dat.motif .= string(i) #label motif
    motif_sim_dat = invasive_sim(inv_motif_dat,missing,"auto",repeat([0.1],length(inv_motif_dat[i]["Neq"])),(0,300),signal_invasive,100,noise_col = "white",saveat=1.0) #simulate stochastic communities.
    #only difference to above is that noise scalar is set generically to be length of number of species
    motif_csv = extract_ts_for_csv(motif_sim_dat,string(i),100:300,5) #crop to t100 - t400 and label as motif "1"
    CSV.write(string("Data/simulations/15_spp/invasive/stress/motif_",i,".csv"),motif_csv,compress=true) #save simulations
    append!(out_df,motif_jac_dat)

    motif_sim_unstress_dat = invasive_sim(inv_motif_dat,missing,"auto",repeat([0.1],length(inv_motif_dat[i]["Neq"])),(0,300),signal_null,100,noise_col = "white",saveat=1.0) #simulate stochastic communities.
    motif_unstress_csv = extract_ts_for_csv(motif_sim_unstress_dat,string(i),100:300,5) #crop to t100 - t400 and label as motif "1"
    CSV.write(string("Data/simulations/15_spp/invasive/unstressed/motif_",i,".csv"),motif_unstress_csv,compress=true) #save simulations

if i == last(motifs)   
    CSV.write(string("Data/simulations/15_spp/invasive/stress/jacobian_15_spp.csv"),out_df,compress=false) #save out jacobian info on last iteration
end
end

################################################
## 20 Species ##
################################################
#Harvesting model auto
    out_df =DataFrame([[],[],[],[],[]], [:community,:collapse,:stress_param,:max_eigval,:motif]) #array to store jacobian information

    for i in motifs
        #i = 1
        motif_dat = RData.load(string("Data/networks/20_spp/web_",i,".RData"),convert=true)["out_file"] #load motif networks generated in R
        #motif_jac_dat = extract_max_eigvec(motif_dat,0.0:0.01:3.0, 100, "auto"; model = "harvest",abs_max = false,concrete_jac = true) #simulate deterministic model for each community
        motif_jac_dat = extract_max_eigvec_alt(motif_dat[1:100],0.0:0.01:ctrl_harvest, 100, "auto"; model = "harvest",abs_max = false) #simulate deterministic model for each community
        motif_jac_dat.motif .= string(i) #label motif
        motif_sim_dat = motif_sim(motif_dat[1:100],missing,"auto",repeat([0.1],length(motif_dat[i]["Neq"])),(0,300),signal_harvest,100,noise_col = "white",saveat=1.0)  #simulate stochastic communities
        #in order, arguments are: parameter data (interaction matrix, starting abundances, carrying capacities etc), starting abundance (if "missing", uses data from the previous argument),
        #species to stress (if "auto", targets the most abundant species in the third trophic level), noise scalar, start-endpoint of timeseries,
        #stress vector, number of simulations (per community), noise colour (white or correlated), saveat what resolution (dt = 0.1) 
        # motif_sim_dat = motif_sim(motif_dat,harvested_spp_ls[string(i)],[0.1,0.1,0.1,0.1],(0,500),signal_motif,100,noise_col = "white",saveat=1.0) 
        motif_csv = extract_ts_for_csv(motif_sim_dat,string(i),100:300,5) #crop to t100 - t300, label as motif "i" and round to 5 decimal places
        CSV.write(string("Data/simulations/20_spp/harvest/stress/motif_",i,".csv"),motif_csv,compress=true) #save simulations
        append!(out_df,motif_jac_dat) #add jacobian info to expanding dataframe

        motif_sim_unstress_dat = motif_sim(motif_dat[1:100],missing,"auto",repeat([0.1],length(motif_dat[i]["Neq"])),(0,300),signal_null,100,noise_col = "white",saveat=1.0)  #simulate stochastic communities but unstressed
        motif_unstress_csv = extract_ts_for_csv(motif_sim_unstress_dat,string(i),100:300,5) 
        CSV.write(string("Data/simulations/20_spp/harvest/unstressed/motif_",i,".csv"),motif_unstress_csv,compress=true) #save simulations

    if i == last(motifs)   
        CSV.write(string("Data/simulations/20_spp/harvest/stress/jacobian_20_spp.csv"),out_df,compress=false) #save out jacobian info on last iteration
    end
    end

#Harvesting model 20%
out_df =DataFrame([[],[],[],[],[]], [:community,:collapse,:stress_param,:max_eigval,:motif]) #array to store jacobian information
harv_spp = [17,18,19,20]

    for i in motifs
        #i = 1
        motif_dat = RData.load(string("Data/networks/20_spp/web_",i,".RData"),convert=true)["out_file"] #load motif networks generated in R
        #motif_jac_dat = extract_max_eigvec(motif_dat,0.0:0.01:3.0, 100, "auto"; model = "harvest",abs_max = false,concrete_jac = true) #simulate deterministic model for each community
        motif_jac_dat = extract_max_eigvec_alt(motif_dat[1:100],0.0:0.01:ctrl_harvest, 100, harv_spp; model = "harvest",abs_max = false) #simulate deterministic model for each community
        motif_jac_dat.motif .= string(i) #label motif
        motif_sim_dat = motif_sim(motif_dat[1:100],missing,harv_spp,repeat([0.1],length(motif_dat[i]["Neq"])),(0,300),signal_harvest,100,noise_col = "white",saveat=1.0)  #simulate stochastic communities
        #in order, arguments are: parameter data (interaction matrix, starting abundances, carrying capacities etc), starting abundance (if "missing", uses data from the previous argument),
        #species to stress (if "auto", targets the most abundant species in the third trophic level), noise scalar, start-endpoint of timeseries,
        #stress vector, number of simulations (per community), noise colour (white or correlated), saveat what resolution (dt = 0.1) 
        # motif_sim_dat = motif_sim(motif_dat,harvested_spp_ls[string(i)],[0.1,0.1,0.1,0.1],(0,500),signal_motif,100,noise_col = "white",saveat=1.0) 
        motif_csv = extract_ts_for_csv(motif_sim_dat,string(i),100:300,5) #crop to t100 - t300, label as motif "i" and round to 5 decimal places
        CSV.write(string("Data/simulations/20_spp/harvest20/stress/motif_",i,"_20percent.csv"),motif_csv,compress=true) #save simulations
        append!(out_df,motif_jac_dat) #add jacobian info to expanding dataframe

        motif_sim_unstress_dat = motif_sim(motif_dat[1:100],missing,"auto",repeat([0.1],length(motif_dat[i]["Neq"])),(0,300),signal_null,100,noise_col = "white",saveat=1.0)  #simulate stochastic communities but unstressed
        motif_unstress_csv = extract_ts_for_csv(motif_sim_unstress_dat,string(i),100:300,5) 
        CSV.write(string("Data/simulations/20_spp/harvest20/unstressed/motif_",i,"_20percent.csv"),motif_unstress_csv,compress=true) #save simulations

    if i == last(motifs)   
        CSV.write(string("Data/simulations/20_spp/harvest20/stress/jacobian_20_spp_20percent.csv"),out_df,compress=false) #save out jacobian info on last iteration
    end
    end

    #Invasive
    out_df =DataFrame([[],[],[],[],[]], [:community,:collapse,:stress_param,:max_eigval,:motif]) #array to store jacobian information

    for i in motifs
        motif_dat = RData.load(string("Data/networks/20_spp/web_",i,".RData"),convert=true)["out_file"] #load motif networks generated in R

        inv_motif_dat = Vector{DictoVec{}}(undef,length(motif_dat[1:100])) #invasive model requires extra community info. Prepare array to hold this info
    for j in eachindex(motif_dat[1:100])

        inv_motif_dat[j] = RData.DictoVec([-1*motif_dat[j]["A_matrix"],
        rand(Distributions.Uniform(0.9,1.1),length(motif_dat[j]["Neq"])),
        rand(Distributions.Uniform(5,15),length(motif_dat[j]["Neq"])),
        repeat([0.1],length(motif_dat[j]["Neq"])),
        #rand(Distributions.Uniform(0,1),length(motif_dat[j]["Neq"])), #random eta
        repeat([1.0],length(motif_dat[j]["Neq"])), #set eta
        motif_dat[j]["Neq"]*10,
        motif_dat[j]["tlvl"]],
        ["A_matrix","R","K","imm","eta","Neq","tlvl"]) #generate additional random info. A_matrix = loaded interaction matrix, R = growth rate,
        #K = carrying capacity, imm = immigration term, eta = stress scalar, Neq = starting abundance, tlvl = trophic level of each species
    end

        motif_jac_dat = extract_max_eigvec_alt(inv_motif_dat,0.0:0.01:ctrl_invasive, 100, "auto"; model = "invasive",abs_max = false) #simulate deterministic model for each community
        motif_jac_dat.motif .= string(i) #label motif
        motif_sim_dat = invasive_sim(inv_motif_dat,missing,"auto",repeat([0.1],length(inv_motif_dat[i]["Neq"])),(0,300),signal_invasive,100,noise_col = "white",saveat=1.0) #simulate stochastic communities.
        #only difference to above is that noise scalar is set generically to be length of number of species
        motif_csv = extract_ts_for_csv(motif_sim_dat,string(i),100:300,5) #crop to t100 - t400 and label as motif "1"
        CSV.write(string("Data/simulations/20_spp/invasive/stress/motif_",i,".csv"),motif_csv,compress=true) #save simulations
        append!(out_df,motif_jac_dat)

        motif_sim_unstress_dat = invasive_sim(inv_motif_dat,missing,"auto",repeat([0.1],length(inv_motif_dat[i]["Neq"])),(0,300),signal_null,100,noise_col = "white",saveat=1.0) #simulate stochastic communities.
        motif_unstress_csv = extract_ts_for_csv(motif_sim_unstress_dat,string(i),100:300,5) #crop to t100 - t400 and label as motif "1"
        CSV.write(string("Data/simulations/20_spp/invasive/unstressed/motif_",i,".csv"),motif_unstress_csv,compress=true) #save simulations

    if i == last(motifs)   
        CSV.write(string("Data/simulations/20_spp/invasive/stress/jacobian_20_spp.csv"),out_df,compress=false) #save out jacobian info on last iteration
    end
    end

################################################
## 25 Species ##
################################################
#Harvesting model auto
out_df =DataFrame([[],[],[],[],[]], [:community,:collapse,:stress_param,:max_eigval,:motif]) #array to store jacobian information

for i in motifs
    #i = 1
    motif_dat = RData.load(string("Data/networks/25_spp/web_",i,".RData"),convert=true)["out_file"] #load motif networks generated in R
    #motif_jac_dat = extract_max_eigvec(motif_dat,0.0:0.01:3.0, 100, "auto"; model = "harvest",abs_max = false,concrete_jac = true) #simulate deterministic model for each community
    motif_jac_dat = extract_max_eigvec_alt(motif_dat[1:100],0.0:0.01:ctrl_harvest, 100, "auto"; model = "harvest",abs_max = false) #simulate deterministic model for each community
    motif_jac_dat.motif .= string(i) #label motif
    motif_sim_dat = motif_sim(motif_dat[1:100],missing,"auto",repeat([0.1],length(motif_dat[i]["Neq"])),(0,300),signal_harvest,100,noise_col = "white",saveat=1.0)  #simulate stochastic communities
    #in order, arguments are: parameter data (interaction matrix, starting abundances, carrying capacities etc), starting abundance (if "missing", uses data from the previous argument),
    #species to stress (if "auto", targets the most abundant species in the third trophic level), noise scalar, start-endpoint of timeseries,
    #stress vector, number of simulations (per community), noise colour (white or correlated), saveat what resolution (dt = 0.1) 
    # motif_sim_dat = motif_sim(motif_dat,harvested_spp_ls[string(i)],[0.1,0.1,0.1,0.1],(0,500),signal_motif,100,noise_col = "white",saveat=1.0) 
    motif_csv = extract_ts_for_csv(motif_sim_dat,string(i),100:300,5) #crop to t100 - t300, label as motif "i" and round to 5 decimal places
    CSV.write(string("Data/simulations/25_spp/harvest/stress/motif_",i,".csv"),motif_csv,compress=true) #save simulations
    append!(out_df,motif_jac_dat) #add jacobian info to expanding dataframe

    motif_sim_unstress_dat = motif_sim(motif_dat[1:100],missing,"auto",repeat([0.1],length(motif_dat[i]["Neq"])),(0,300),signal_null,100,noise_col = "white",saveat=1.0)  #simulate stochastic communities but unstressed
    motif_unstress_csv = extract_ts_for_csv(motif_sim_unstress_dat,string(i),100:300,5) 
    CSV.write(string("Data/simulations/25_spp/harvest/unstressed/motif_",i,".csv"),motif_unstress_csv,compress=true) #save simulations

if i == last(motifs)   
    CSV.write(string("Data/simulations/25_spp/harvest/stress/jacobian_25_spp.csv"),out_df,compress=false) #save out jacobian info on last iteration
end
end

#Harvesting model 20%
out_df =DataFrame([[],[],[],[],[]], [:community,:collapse,:stress_param,:max_eigval,:motif]) #array to store jacobian information
harv_spp = [21,22,23,24,25]
for i in motifs
    #i = 1
    motif_dat = RData.load(string("Data/networks/25_spp/web_",i,".RData"),convert=true)["out_file"] #load motif networks generated in R
    #motif_jac_dat = extract_max_eigvec(motif_dat,0.0:0.01:3.0, 100, "auto"; model = "harvest",abs_max = false,concrete_jac = true) #simulate deterministic model for each community
    motif_jac_dat = extract_max_eigvec_alt(motif_dat[1:100],0.0:0.01:ctrl_harvest, 100, harv_spp; model = "harvest",abs_max = false) #simulate deterministic model for each community
    motif_jac_dat.motif .= string(i) #label motif
    motif_sim_dat = motif_sim(motif_dat[1:100],missing,harv_spp,repeat([0.1],length(motif_dat[i]["Neq"])),(0,300),signal_harvest,100,noise_col = "white",saveat=1.0)  #simulate stochastic communities
    #in order, arguments are: parameter data (interaction matrix, starting abundances, carrying capacities etc), starting abundance (if "missing", uses data from the previous argument),
    #species to stress (if "auto", targets the most abundant species in the third trophic level), noise scalar, start-endpoint of timeseries,
    #stress vector, number of simulations (per community), noise colour (white or correlated), saveat what resolution (dt = 0.1) 
    # motif_sim_dat = motif_sim(motif_dat,harvested_spp_ls[string(i)],[0.1,0.1,0.1,0.1],(0,500),signal_motif,100,noise_col = "white",saveat=1.0) 
    motif_csv = extract_ts_for_csv(motif_sim_dat,string(i),100:300,5) #crop to t100 - t300, label as motif "i" and round to 5 decimal places
    CSV.write(string("Data/simulations/25_spp/harvest20/stress/motif_",i,"_20percent.csv"),motif_csv,compress=true) #save simulations
    append!(out_df,motif_jac_dat) #add jacobian info to expanding dataframe

    motif_sim_unstress_dat = motif_sim(motif_dat[1:100],missing,"auto",repeat([0.1],length(motif_dat[i]["Neq"])),(0,300),signal_null,100,noise_col = "white",saveat=1.0)  #simulate stochastic communities but unstressed
    motif_unstress_csv = extract_ts_for_csv(motif_sim_unstress_dat,string(i),100:300,5) 
    CSV.write(string("Data/simulations/25_spp/harvest20/unstressed/motif_",i,"_20percent.csv"),motif_unstress_csv,compress=true) #save simulations

if i == last(motifs)   
    CSV.write(string("Data/simulations/25_spp/harvest20/stress/jacobian_25_spp_20percent.csv"),out_df,compress=false) #save out jacobian info on last iteration
end
end

#Invasive
out_df =DataFrame([[],[],[],[],[]], [:community,:collapse,:stress_param,:max_eigval,:motif]) #array to store jacobian information

for i in motifs
    motif_dat = RData.load(string("Data/networks/25_spp/web_",i,".RData"),convert=true)["out_file"] #load motif networks generated in R

    inv_motif_dat = Vector{DictoVec{}}(undef,length(motif_dat[1:100])) #invasive model requires extra community info. Prepare array to hold this info
for j in eachindex(motif_dat[1:100])

    inv_motif_dat[j] = RData.DictoVec([-1*motif_dat[j]["A_matrix"],
    rand(Distributions.Uniform(0.9,1.1),length(motif_dat[j]["Neq"])),
    rand(Distributions.Uniform(5,15),length(motif_dat[j]["Neq"])),
    repeat([0.1],length(motif_dat[j]["Neq"])),
    #rand(Distributions.Uniform(0,1),length(motif_dat[j]["Neq"])), #random eta
    repeat([1.0],length(motif_dat[j]["Neq"])), #set eta
    motif_dat[j]["Neq"]*10,
    motif_dat[j]["tlvl"]],
    ["A_matrix","R","K","imm","eta","Neq","tlvl"]) #generate additional random info. A_matrix = loaded interaction matrix, R = growth rate,
    #K = carrying capacity, imm = immigration term, eta = stress scalar, Neq = starting abundance, tlvl = trophic level of each species
end

    motif_jac_dat = extract_max_eigvec_alt(inv_motif_dat,0.0:0.01:ctrl_invasive, 100, "auto"; model = "invasive",abs_max = false) #simulate deterministic model for each community
    motif_jac_dat.motif .= string(i) #label motif
    motif_sim_dat = invasive_sim(inv_motif_dat,missing,"auto",repeat([0.1],length(inv_motif_dat[i]["Neq"])),(0,300),signal_invasive,100,noise_col = "white",saveat=1.0) #simulate stochastic communities.
    #only difference to above is that noise scalar is set generically to be length of number of species
    motif_csv = extract_ts_for_csv(motif_sim_dat,string(i),100:300,5) #crop to t100 - t400 and label as motif "1"
    CSV.write(string("Data/simulations/25_spp/invasive/stress/motif_",i,".csv"),motif_csv,compress=true) #save simulations
    append!(out_df,motif_jac_dat)

    motif_sim_unstress_dat = invasive_sim(inv_motif_dat,missing,"auto",repeat([0.1],length(inv_motif_dat[i]["Neq"])),(0,300),signal_null,100,noise_col = "white",saveat=1.0) #simulate stochastic communities.
    motif_unstress_csv = extract_ts_for_csv(motif_sim_unstress_dat,string(i),100:300,5) #crop to t100 - t400 and label as motif "1"
    CSV.write(string("Data/simulations/25_spp/invasive/unstressed/motif_",i,".csv"),motif_unstress_csv,compress=true) #save simulations

if i == last(motifs)   
    CSV.write(string("Data/simulations/25_spp/invasive/stress/jacobian_25_spp.csv"),out_df,compress=false) #save out jacobian info on last iteration
end
end
