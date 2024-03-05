using ArgParse
using GerryChain
using Statistics
using ResumableFunctions


# From https://github.com/mggg/GerryChainJulia/pull/123
@resumable function short_bursts_recom_iter(
    score::AbstractScore,
    burst_length::Int,
    num_bursts::Int,
    graph::BaseGraph,
    partition::Partition,
    pop_constraint::PopulationConstraint,
    scores::Array{S,1};
    acceptance_fn::F = always_accept,
    verbose::Bool = false,
) where {F<:Function, S<:AbstractScore}
    best_partition = partition
    best_score = deepcopy(eval_score_on_partition(graph, partition, score))
    if verbose
        iter = ProgressBar(1:num_bursts)
    else
        iter = 1:num_bursts
    end

    iter_function = recom_chain_iter
    for _ = iter
        for (partition, score_vals) in iter_function(
            graph,
            deepcopy(best_partition),
            pop_constraint,
            burst_length,
            cat(score, scores, dims=1),
            acceptance_fn = acceptance_fn,
            progress_bar = false
        )
            if score_vals[score.name] >= best_score
                best_partition = deepcopy(partition)
                best_score = deepcopy(eval_score_on_partition(graph, partition, score))
            end
            @yield partition, score_vals
        end
    end
end

s = ArgParseSettings()
@add_arg_table s begin
	"state"
		help = "Abbreviation of state (i.e. RI)"
		arg_type = String
		required = true
	"num_steps"
		help = "Number of steps to run Markov chain"
		arg_type = Int
		required = true
	"accept_fn"
		help = "Which acceptance function to use when generating ensemble"
		arg_type = String
		range_tester = (x->x=="all"||x=="burst")
		required = true
	"pop_dt"
		help = "Type of population data (dp, sf)"
		arg_type = String
		range_tester = (x->x=="dp"||x=="sf")
		required = true
	"level"
		help = "Level of geography"
		arg_type = String
		range_tester = (x->x=="block"||x=="bg"||x=="tract")
		required = true
	"district_type"
		help = "Which districts to test (State House, State Senate, Congressional, Block Group)"
		arg_type = String
		range_tester = (x->x=="SH"||x=="SS"||x=="CD"||x=="BG")
		required = true
	"group"
		help = "Which group to optimize (white, black, hisp)"
		arg_type = String
		range_tester = (x->x=="white"||x=="black"||x=="hisp")
		required = true
    "exp_name"
        help = "Name of experiment (will create directory with this name to store results)"
        arg_type = String
        required = true
    "--shp_path"
        help = "Path to directory containing preprocessed shapefiles"
        arg_type = String
        default = "./data/"
        required = false
    "--save_dir"
        help = "Path to save ensembles"
        arg_type = String
        default = "./chains/"
        required = false
	"--threshold"
        help = "Threshold for effective district proportion"
        arg_type = Float64
        range_tester = (x->x>=0.0&&x<=1.0)
        default = 0.5
        required = false
	"--burst_length"
        help = "Length of bursts"
        arg_type = Int
        range_tester = (x->x>0)
        default = 10
        required = false
    "--num_chains"
        help="Number of chains to run"
        arg_type = Int
        range_tester = (x->x>0)
        default = 10
        required = false
    "--data_vintage"
        help = "Vintage of demo data to use"
        arg_type = String
        default = "20210608"
        range_tester = (x->x=="20210608"||x=="20200527")
        required = false
    "--epsilon"
        help = "Allowed deviation from desired district population"
        arg_type = Float64
        default = 0.1
        range_tester = (x->x>=0.0)
        required = false
    "--totpop"
        help = "Log total district populations (sf and dp)"
        action = :store_true
    "--thresh_name"
        help = "Include threshold in output file name"
        action = :store_true
    "--nowrite"
        help = "Don't write output to file"
        action = :store_true
    "--nosave"
        help = "Don't save score csv"
        action = :store_true
    "--burn_time"
        help = "Number of steps to burn in"
        arg_type = Int
        default = 0
        range_tester = (x->x>=0)
        required = false
    "--vap"
        help = "Use VAP numbers when computing majority-minority districts"
        action = :store_true
    "--save_every"
        help = "Rate at which to subsample data (i.e. write every nth plan)"
        arg_type = Int
        default = 1
    "--detailed"
        help = "Save detailed demographic data"
        action = :store_true
    "--suffix"
        help = "Suffix for output file name"
        arg_type = String
        default = ""
end

args = parse_args(ARGS, s)
state = args["state"]
num_steps = args["num_steps"]
accept_fn = args["accept_fn"]
level = args["level"]
district_type = args["district_type"]
group = args["group"]
effective_threshold = args["threshold"]
burst_length = args["burst_length"]
num_chains = args["num_chains"]
nosave = args["nosave"]
data_vintage = args["data_vintage"]
epsilon = args["epsilon"]
nowrite = args["nowrite"]
exp_name = args["exp_name"]
pop_dt = args["pop_dt"]
totpop = args["totpop"]
thresh_name = args["thresh_name"]
burn_time = args["burn_time"]
vap = args["vap"]
save_every = args["save_every"]
detailed = args["detailed"]
suffix = args["suffix"]


# pop_dt = "dp"
if pop_dt == "sf"
    noise_type = "dp"
else
    noise_type = pop_dt
end

base_path = args["shp_path"]
save_dir = args["save_dir"]

graph_path = base_path*"merged_shp_$(data_vintage)/$(state)/$(level)/$(state)_$(level)_$(district_type).shp"

pop_col = "TOTPOP"*pop_dt
assignment_col = "DISTRICT"

# Initialize graph and partition
graph = BaseGraph(graph_path, pop_col)
partition = Partition(graph, assignment_col)

# Define parameters of chain (number of steps and population constraint)
pop_constraint = PopulationConstraint(graph, partition, epsilon)

# Initialize Election of interest
if vap 
    if group == "white"
        ename = "PROPWHITE"
        name1 = "WVAP"
        name2 = "NWVAP"
        short = "w"

        score_name = "nw"
        score_col = "NWVAP"

        nscore_name = "w"
        nscore_col = "WVAP"
    elseif group == "black"
        ename = "PROPBLACK"
        name1 = "BVAP"
        name2 = "NBVAP"
        short = "b"

        score_name = "b"
        score_col = "BVAP"
        nscore_name = "nb"
        nscore_col = "NBVAP"
    elseif group == "hisp"
        ename = "PROPHISP"
        name1 = "HVAP"
        name2 = "NHVAP"
        short = "h"

        score_name = "h"
        score_col = "HVAP"
        nscore_name = "nh"
        nscore_col = "NHVAP"
    else
        error("Invalid group")
    end
else
    if group == "white"
        ename = "PROPWHITE"
        name1 = "NH_WHITE"
        name2 = "NWHITE"
        short = "w"

        score_name = "nw"
        score_col = "NWHITE"
        nscore_name = "w"
        nscore_col = "NH_WHITE"
    elseif group == "black"
        ename = "PROPBLACK"
        name1 = "NH_BLACK"
        name2 = "NBLACK"
        short = "b"

        score_name = "b"
        score_col = "NH_BLACK"
        nscore_name = "nb"
        nscore_col = "NBLACK"
    elseif group == "hisp"
        ename = "PROPHISP"
        name1 = "HISP"
        name2 = "NHISP"
        short = "h"

        score_name = "h"
        score_col = "HISP"
        nscore_name = "nh"
        nscore_col = "NHISP"
    else
        error("Invalid group")
    end
end

electiondp = Election(ename*noise_type, [name1*noise_type, name2*noise_type], partition.num_dists)
electionsf = Election(ename*"sf", [name1*"sf", name2*"sf"], partition.num_dists)

if pop_dt == noise_type
	election = electiondp
else
	election = electionsf
end

# Define election-related metrics and scores
election_metrics = [
    vote_count("count_$(short)"*noise_type, electiondp, name1*noise_type),
    vote_count("count_$(short)_sf", electionsf, name1*"sf"),
]

# Score plans based on number of majority-minority districts
function get_seats_score(dt::String, sname::String, scol::String, non_sname::String, non_scol::String, name::String)
	function score_fn(graph::BaseGraph, partition::Partition)
		group = DistrictAggregate(sname, scol*dt)
		num_group = [eval_score_on_district(graph, partition, group, i) for i=1:partition.num_dists]
		nongroup = DistrictAggregate(non_sname, non_scol*dt)
		num_nongroup = [eval_score_on_district(graph, partition, nongroup, i) for i=1:partition.num_dists]
        return sum([num_group[i] / (num_group[i] + num_nongroup[i]) > effective_threshold for i=1:partition.num_dists if num_group[i] + num_nongroup[i] > 0])
	end
	return PlanScore(name, score_fn)
end

scores = [
    get_seats_score(noise_type, score_name*"_"*noise_type, score_col, nscore_name*"_"*noise_type, nscore_col, "num_$(score_name)_seats_"*noise_type),
	get_seats_score("sf", score_name*"_sf", score_col, nscore_name*"_sf", nscore_col, "num_$(score_name)_seats_sf")
]

if vap 
    if detailed 
        scores = cat(scores, Array([DistrictAggregate("VAPsf"), DistrictAggregate("VAP"*noise_type), 
                                    DistrictAggregate("WVAPsf"), DistrictAggregate("WVAP"*noise_type),
                                    DistrictAggregate("BVAPsf"), DistrictAggregate("BVAP"*noise_type), 
                                    DistrictAggregate("HVAPsf"), DistrictAggregate("HVAP"*noise_type),
                                    DistrictAggregate("ASIANVAPsf"), DistrictAggregate("ASIANVAP"*noise_type), 
                                    DistrictAggregate("AMINVAPsf"), DistrictAggregate("AMINVAP"*noise_type),
                                    DistrictAggregate("NHPIVAPsf"), DistrictAggregate("NHPIVAP"*noise_type)]
                                    ), dims=1)
    else
        scores = cat(scores, Array([DistrictAggregate(name2*"sf"), DistrictAggregate(name2*noise_type), DistrictAggregate(name1*"sf"), DistrictAggregate(name1*noise_type)]), dims=1)
    end
else
    if detailed 
        scores = cat(scores, Array([DistrictAggregate("TOTPOPsf"), DistrictAggregate("TOTPOP"*noise_type), 
                                    DistrictAggregate("NH_WHITEsf"), DistrictAggregate("NH_WHITE"*noise_type),
                                    DistrictAggregate("NH_BLACKsf"), DistrictAggregate("NH_BLACK"*noise_type), 
                                    DistrictAggregate("HISPsf"), DistrictAggregate("HISP"*noise_type),
                                    DistrictAggregate("NH_ASIANsf"), DistrictAggregate("NH_ASIAN"*noise_type), 
                                    DistrictAggregate("NH_AMINsf"), DistrictAggregate("NH_AMIN"*noise_type),
                                    DistrictAggregate("NH_NHPIsf"), DistrictAggregate("NH_NHPI"*noise_type)]
                                    ), dims=1)
    else
        scores = cat(scores, Array([DistrictAggregate("TOTPOPsf"), DistrictAggregate("TOTPOP"*noise_type), DistrictAggregate(name1*"sf"), DistrictAggregate(name1*noise_type)]), dims=1)
    end
end
if totpop
    scores = cat(scores, Array([DistrictAggregate("TOTPOPsf"), DistrictAggregate("TOTPOP"*noise_type)]), dims=1)
end

# Define scoring function for short bursts
election_metrics = [
    vote_count("count_$(score_name)_"*noise_type, electiondp, score_col*noise_type),
    vote_share("share_$(score_name)_"*noise_type, electiondp, score_col*noise_type),
    vote_count("count_$(score_name)_sf", electionsf, score_col*"sf"),
    vote_share("share_$(score_name)_sf", electionsf, score_col*"sf"),
]
election_score = ElectionTracker(election, election_metrics)
function score_group(graph::BaseGraph, partition::Partition)
    district_scores = eval_score_on_partition(graph, partition, election_score)
    group_share = district_scores["share_$(score_name)_$(pop_dt)"]
    minority_dists = filter(x -> x < effective_threshold, group_share)
    if isempty(minority_dists)
        reward = 0
    else
        reward = maximum(minority_dists)
    end
    score = length(filter(x -> x >= effective_threshold, group_share)) + reward
    return score
end
reward_partial = PlanScore("reward_$(score_name)", score_group)

if burn_time > 0
    if nowrite
        println("Running burn-in")
    end
    for (partition, score_vals) in recom_chain_iter(graph, partition, pop_constraint, burn_time, scores, progress_bar=false)
        continue
    end
end
original_partition = deepcopy(partition)

first_scores = score_initial_partition(graph, original_partition, scores)
chain_scores = ChainScoreData(deepcopy(scores), [first_scores])

if nowrite
    println("Running chains")
end
if accept_fn == "all"
    chain_function = recom_chain_iter
    for (plan_num, (part, score_vals)) in enumerate(chain_function(graph, partition, pop_constraint, num_steps, scores, progress_bar=false))
        if !nowrite
            write(stdout, repr(part.assignments .- 1), "\n")
        end
        push!(chain_scores.step_values, score_vals)
    end
else
    steps_per_chain = num_steps ÷ num_chains
    num_bursts = steps_per_chain ÷ burst_length
    for chain_num in 1:num_chains
        orig_part_copy = deepcopy(original_partition)
        for (plan_num, (part, score_vals)) in enumerate(short_bursts_recom_iter(reward_partial, burst_length, num_bursts, graph, orig_part_copy, pop_constraint, scores))
            if !nowrite
                write(stdout, repr(part.assignments .- 1), "\n")
            end

            push!(chain_scores.step_values, score_vals)
        end
    end
end

flush(stdout)

loc_name = "$(state)"

if detailed 
    detailed_str = "_detailed"
else
    detailed_str = ""
end

if !nosave
    if accept_fn == "all"
        save_dir = save_dir*"chain_data/$(exp_name)/$(data_vintage)/"
        mkpath(save_dir)
        save_name = save_dir * "$(loc_name)_$(district_type)_$(level)_$(num_steps)_$(pop_dt)_all$(detailed_str)$(suffix).csv"
        if thresh_name
            save_name = save_dir * "$(loc_name)_$(district_type)_$(level)_$(num_steps)_$(pop_dt)_all$(detailed_str)_$(epsilon)$(suffix).csv"
        end
        save_scores_to_csv(save_name, chain_scores, String[], save_every)
    elseif accept_fn == "burst"
        num_bursts = num_steps ÷ burst_length
        save_dir = save_dir*"short_burst_chain_data/$(exp_name)/$(data_vintage)/"
        mkpath(save_dir)
        save_name = save_dir * "$(loc_name)_$(district_type)_$(level)_$(burst_length)_$(num_bursts)_$(pop_dt)_bursts$(detailed_str)$(suffix).csv"
        if thresh_name
            save_name = save_dir * "$(loc_name)_$(district_type)_$(level)_$(burst_length)_$(num_bursts)_$(pop_dt)_bursts$(detailed_str)_$(epsilon)$(suffix).csv"
        end
        save_scores_to_csv(save_name, chain_scores, String[], save_every)
    end
end
