# Main program for the 5G NR Simulation

include("Simulation_5G.jl")
include("Protocols_5G.jl")
using MAT
max_packet_rate = 256
λ_arr = 1:200

source = Source(1,max_packet_rate,100)
# It is a traffic shaper that can handle bursty traffic
source_out = Performance(source)[1]

sdpa_l = SDAP(source_out)
sdap_out = Performance(sdpa_l)[1]

pdcp_l = PDCP(sdap_out)
pdcp_out = Performance(pdcp_l)[1]

rlc_l = RLC_AM(pdcp_out)
rlc_out = Performance(rlc_l)[1]

mac_l = MAC(rlc_out)
mac_out = Performance(mac_l)[1]

phy_l = PHY(mac_out)

stack = CrossLayer(source,sdpa_l,pdcp_l,rlc_l,mac_l,phy_l)

# We want the event number will be bounded
event_no = 500000

# The final simulation time should be derived from this value

TTI = 0.0005
for λ in λ_arr
    sim_time_limit = event_no/λ
    stack.S.λ_raw = λ
    cls = CrossLayer_Simulation(nothing,nothing,0.0,sim_time_limit,[],[],(1/λ),stack,[],false,Inf,[],false,Inf,[],false,Inf,[],[],Inf,[],[],false,Inf,[],TTI)
    print("Current Arrival Rate is: ")
    println(λ)
    Simulation_Cycle(cls)
end
