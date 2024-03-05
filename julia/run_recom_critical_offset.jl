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
    # end
end

s = ArgParseSettings()
@add_arg_table s begin
	"state"
		help = "Abbreviation of state (i.e. RI)"
		arg_type = String
		required = true
	"num_steps"
		help = "Number of plans in ensemble"
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
		help = "Which districts to test (State House, State Senate, Congressional)"
		arg_type = String
		range_tester = (x->x=="SH"||x=="SS"||x=="CD"||x=="BG")
		required = true
	"group"
		help = "Which group to optimize (white, black, hisp)"
		arg_type = String
		range_tester = (x->x=="white"||x=="black"||x=="hisp")
		required = true
    "--shp_path"
        help = "Path to directory containing preprocessed shapefiles"
        arg_type = String
        default = "./data/"
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
    "--thresh_name"
        help = "Include threshold in output file name"
        action = :store_true
    "--nowrite"
        help = "Don't write output to file"
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
    "--detailed"
        help = "Save detailed demographic data"
        action = :store_true
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
data_vintage = args["data_vintage"]
nowrite = args["nowrite"]
pop_dt = args["pop_dt"]
thresh_name = args["thresh_name"]
burn_time = args["burn_time"]
vap = args["vap"]
detailed = args["detailed"]

# pop_dt = "dp"
if pop_dt == "sf"
    noise_type = "dp"
else
    noise_type = pop_dt
end

base_path = args["shp_path"]

graph_path = base_path*"merged_shp_$(data_vintage)/$(state)/$(level)/$(state)_$(level)_$(district_type).shp"

pop_col = "TOTPOP"*pop_dt
assignment_col = "DISTRICT"

# Initialize graph and partition
graph = BaseGraph(graph_path, pop_col)
partition = Partition(graph, assignment_col)

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

scores = Array([DistrictAggregate("TOTPOPsf"), DistrictAggregate("TOTPOP"*noise_type)])

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

println("pop_tol\tDiscrep. Rate")
for epsilon in [0.05, 0.0495, 0.049, 0.0485, 0.048, 0.0475, 0.047, 0.046, 0.0455, 0.045, 0.0445, 0.044, 0.0435, 0.043, 0.0425, 0.042, 0.0415, 0.041, 0.0405, 0.04]
    # Define parameters of chain (number of steps and population constraint)
    pop_constraint = PopulationConstraint(graph, partition, epsilon)

    if burn_time > 0
        for (partition, score_vals) in recom_chain_iter(graph, partition, pop_constraint, burn_time, scores, progress_bar=false)
            continue
        end
    end
    original_partition = deepcopy(partition)

    first_scores = score_initial_partition(graph, original_partition, scores)
    chain_scores = ChainScoreData(deepcopy(scores), [first_scores])

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

    # find proportion of plans that exceed pop constraint in SF data
    # break if proportion is less then 2%
    pop_data = get_score_values(chain_scores, "TOTPOPsf")
    pop_max = maximum(pop_data, dims=2)
    pop_min = minimum(pop_data, dims=2)
    pop_mean = mean(pop_data, dims=2)
    pop_discrep = (pop_max - pop_min) ./ pop_mean
    prop_exceed = sum(pop_discrep .> 0.1) / num_steps
    println("$(epsilon)\t$(prop_exceed)")
    if prop_exceed < 0.02
        critical_offset = round(0.05 - epsilon; digits=4)
        println("Critical offset: $(critical_offset)")
        break
    end
    if epsilon == 0.04
        println("Critical offset: > 0.01")
    end
end
