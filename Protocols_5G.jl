# Code for Defining Structs for Each of the Sub-layers
using LinearAlgebra

########## SDAP ##########

mutable struct SDAP
    # Struct for defining the parameters for SDAP performance
    λ_SDAP # Arrival Rate to SDAP layer
    μ_SDAP # Service Rate of the SDAP layer
    function SDAP(arrival_rate,service_rate=1/(0.02*(10^-6)))
        return new(arrival_rate,service_rate)   
    end
end

function Performance(L::SDAP)
    # Function for evaluating the performance of SDAP layer
    if L.λ_SDAP < L.μ_SDAP
        Th = min(L.λ_SDAP,L.μ_SDAP)
        Delay = 1/(L.μ_SDAP-L.λ_SDAP)
    else
        Th = L.μ_SDAP
        Delay = 3600.0+(L.λ_SDAP-L.μ_SDAP)
    end
    #if Delay < 0
    #    Delay = Inf
    #end
    return Th,Delay
end

########## PDCP ##########

mutable struct PDCP
    λ_PDPC # Arrival Rate to PDCP layer
    μ1_PDCP # MAC calculing Rate
    μ2_PDCP # Ciphering Rate
    μ3_PDCP # Header Creation Rate
    function PDCP(arrival_rate,MAC_rate=1.140773389596174e+05,Ciph_rate=2.374921261809742e+05,Hdr_rate=1/(0.02*(10^-6)))
        return new(arrival_rate,MAC_rate,Ciph_rate,Hdr_rate)
    end
end


function Performance(L::PDCP)
    if L.λ_PDPC < (L.μ1_PDCP*L.μ2_PDCP*L.μ3_PDCP)/((L.μ2_PDCP*L.μ3_PDCP)+(L.μ1_PDCP*L.μ3_PDCP)+(L.μ1_PDCP*L.μ2_PDCP))
        B00 = [-L.λ_PDPC]
        B01 = [L.λ_PDPC 0 0]
        B10 = [0; 0; L.μ3_PDCP]

        A0 = [0 0 0; 0 0 0; L.μ3_PDCP 0 0]
        A1 = [-(L.λ_PDPC+L.μ1_PDCP) L.μ1_PDCP 0; 0 -(L.λ_PDPC+L.μ2_PDCP) L.μ2_PDCP; 0 0 -(L.λ_PDPC+L.μ3_PDCP)]
        A2 = [L.λ_PDPC 0 0; 0 L.λ_PDPC 0; 0 0 L.λ_PDPC]

        A = A0+A1+A2
        r = 3 # Number of sequential phases
        A[:,r] = ones(r,1)

        rhsA = zeros(1,r)
        rhsA[r] = 1
        
        piA = rhsA/A

        e = ones(r,1)
        mean_arr = piA*A2*e
        mean_dep = piA*A0*e
        if ((mean_arr[1])>=(mean_dep[1]))
            println("Unestable System")
        end

        R = [(L.λ_PDPC*(L.λ_PDPC+L.μ2_PDCP)*(L.λ_PDPC+L.μ3_PDCP))/(L.μ1_PDCP*L.μ2_PDCP*L.μ3_PDCP) (L.λ_PDPC*(L.λ_PDPC+L.μ3_PDCP))/(L.μ2_PDCP*L.μ3_PDCP) L.λ_PDPC/L.μ3_PDCP; ((L.λ_PDPC^2)*(L.λ_PDPC+L.μ2_PDCP+L.μ3_PDCP))/(L.μ1_PDCP*L.μ2_PDCP*L.μ3_PDCP) L.λ_PDPC*(L.λ_PDPC+L.μ3_PDCP)/(L.μ2_PDCP*L.μ3_PDCP) L.λ_PDPC/L.μ3_PDCP; (L.λ_PDPC^2)*(L.λ_PDPC+L.μ2_PDCP)/(L.μ1_PDCP*L.μ2_PDCP*L.μ3_PDCP) (L.λ_PDPC^2)/(L.μ2_PDCP*L.μ3_PDCP) L.λ_PDPC/L.μ3_PDCP]

        Qb = [B00 B01; B10 (A1+(R*A0))]
        Qb[:,r+1] = [1; zeros(r,1)]

        rhs_b = zeros(1,r+1)
        rhs_b[r+1] = 1

        s_limit = rhs_b/Qb

        pi0_nn = [s_limit[1]]
        pi1_nn = transpose(s_limit[2:r+1])

        alpha = pi0_nn[1] + sum(pi1_nn/(Matrix(1.0I,r,r)-R))

        pi0 = pi0_nn/alpha
        pi1 = pi1_nn/alpha

        λ_mean = sum(piA*A2)
        μ_mean = sum(piA*A0)

        EN = sum(pi1/((Matrix(1.0I,r,r)-R)^2))

        ER = EN/λ_mean

        Th = min(λ_mean,μ_mean)
        Delay = ER

        return Th, Delay
    else
        Th_PDCP = (L.μ1_PDCP*L.μ2_PDCP*L.μ3_PDCP)/((L.μ2_PDCP*L.μ3_PDCP)+(L.μ1_PDCP*L.μ3_PDCP)+(L.μ1_PDCP*L.μ2_PDCP))
        Delay_PDCP = 3600.0+(L.μ1_PDCP*L.μ2_PDCP*L.μ3_PDCP)/((L.μ2_PDCP*L.μ3_PDCP)+(L.μ1_PDCP*L.μ3_PDCP)+(L.μ1_PDCP*L.μ2_PDCP))
        
        return Th_PDCP,Delay_PDCP
    end

end

mutable struct RLC_AM
    λ_traffic # Arrival rate
    p # Probability of success in a transmission
    n # Maximum transmission buffer buffer_Size
    max_ret # Maximum number of retransmission
    windows_size # Retransmission window buffer
    λ_header # Header creation rate
    λ_ACK # Successfull transmission rate
    λ_NACK # Erroneous transmission rate
    function RLC_AM(arr_rate,p=0.97,n=20,max_ret=5,windows_size=10,theader=1/(0.02*(10^-6)),ttrx=1/(0.01*(10^-6)))
        return new(arr_rate,p,n,max_ret,windows_size,theader,p*ttrx,(1-p)*ttrx)
    end
end

function Performance(L::RLC_AM)
    RLC_Model = open("RLC_Proposal_Model.txt","r")
    RLC_Model_lines = readlines(RLC_Model)
    close(RLC_Model)
    RLC_Model_lines[5] = "p   "*string(L.p)
    RLC_Model_lines[6] = "n   "*string(L.n)
    RLC_Model_lines[7] = "max_ret   "*string(L.max_ret)
    RLC_Model_lines[8] = "window_size   "*string(L.windows_size)

    RLC_Model_lines[87] = "  Traffic ind "*string(L.λ_traffic)*" guard f_traffic()"
    RLC_Model_lines[88] = "  THeader ind "*string(L.λ_header)
    RLC_Model_lines[89] = "  TACK ind "*string(L.λ_ACK)
    RLC_Model_lines[90] = "  TNACK ind "*string(L.λ_NACK)

    
    string_val = RLC_Model_lines[1]*"\n"
    for ind in 2:length(RLC_Model_lines)
        string_val = string_val*RLC_Model_lines[ind]*"\n"
    end
    RLC_Model = open("RLC_Proposal_Model_modif.txt","w")
    write(RLC_Model,string_val)
    close(RLC_Model)

    # Model evaluation
    command = `C:/Sharpe-Gui/sharpe/sharpe.exe RLC_Proposal_Model_modif.txt`
    results = read(command,String);

    results_split = split(results,"\n")

    buffer_Size_str = results_split[4]

    buffer_Size = parse(Float64,split(buffer_Size_str,"   ")[2])

    Traffic_str = results_split[9]

    Traffic = parse(Float64,split(Traffic_str,"   ")[2])

    TACK_str = results_split[14]

    TACK = parse(Float64,split(TACK_str,"   ")[2])

    SR_Window_str = results_split[19]

    SR_Window = parse(Float64,split(SR_Window_str,"   ")[2])

    THeader_str = results_split[24]

    THeader = parse(Float64,split(THeader_str,"   ")[2])

    Retrans_Buffer_str = results_split[29]

    Retrans_Buffer = parse(Float64,split(Retrans_Buffer_str,"   ")[2])

    TNACK_str = results_split[34]

    TNACK = parse(Float64,split(TNACK_str,"   ")[2])

    Th_RLC_AM = TACK

    Delay_RLC_AM = (buffer_Size/Traffic) + (SR_Window/THeader) +(Retrans_Buffer/TNACK) + L.p/L.λ_ACK


    return Th_RLC_AM, Delay_RLC_AM
end

#################### MAC ####################

mutable struct MAC
    # Struct used to define the MAC layer performance
    λ_MAC # Arrival rate to the MAC layer
    μ_MAC # MAC header creation rate
    TTI # Time Transmission Interval
    channel # Array for Channel
    α # Available Resource 
    μ_TX_MAC # Tranmisión Rate
    function MAC(arrival_rate,service_rate=1/(0.02*(10^-6)),interval=0.5*(10^-3),channel=zeros(14,51),available_rb=0.1)
        channel_size = size(channel)
        tx_mac_rate = (available_rb*(channel_size[1]*channel_size[2]))/interval
        return new(arrival_rate,service_rate,interval,channel,available_rb,tx_mac_rate)
    end
end

function Performance(L::MAC)
    if L.λ_MAC < min(L.μ_MAC,L.μ_TX_MAC)
        Th = min(L.λ_MAC,L.μ_MAC,L.μ_TX_MAC)
        Delay = (1/(L.μ_MAC-L.λ_MAC)) + (L.TTI/2) + (1/(L.μ_TX_MAC-(min(L.λ_MAC,L.μ_MAC))))
    else
        #Th = min(L.μ_MAC,1/L.TTI)
        Th = min(L.λ_MAC,L.μ_MAC,L.μ_TX_MAC)
        Delay = 3600.0+(L.λ_MAC-L.μ_MAC)+(L.TTI/2) + (L.TTI/2)
    end
    return Th,Delay
end

#################### MAC ####################

mutable struct PHY
    λ # Arrival Rate
    J # Number of component carrier
    v_j # Layers for MIMO
    Q_j_m # Modulation Order
    BPSym # Bits per symbol (based on modulation)
    R # Coding Rate
    R_max # Max code rate
    f # Coding factor
    BW_j # Bandwidth
    μ # Numerology
    scs # Sub-carrier spacing (depends on the numerology)
    OH_j # Overhead
    N_BW_j_mu # Number of RB available for the current BW and μ
    T_mu_s # Symbol duration
    N1 #
    N2 #
    Tc #
    k #
    Tproc2
    Tproc1
    radio # Cell radius
    TBS # Transport Block Size
    CQI # Channel Quality Indicator
    function PHY(λ,J=1,v_j=2,Q_j_m=8,R=885,f=1,BW_j=20,μ=1,OH_j=0.08,N_BW_j_mu=51,N1=4.5,N2=5.5,Tc=1/(480000*4096),k=64,radio=866,TBS=2088)
        return new(λ,J,v_j,Q_j_m,2^Q_j_m,R,R/1024,f,BW_j,μ,15*(2^μ),OH_j,N_BW_j_mu,(10^(-3))/(14*(2^μ)),N1,N2,Tc,k,N1*(2048+144)*k*(2.0^-μ)*Tc,max(N2*(2014+144)*k*(2.0^-μ)*Tc,0),radio,TBS,14)
    end
end

function max_PHY_Th(L::PHY)
    Th_bps = 0
    for j in 1:L.J
        Th_bps = Th_bps + L.v_j*L.Q_j_m*L.f*L.R_max*((L.N_BW_j_mu*12)/L.T_mu_s)*(1-L.OH_j)
    end
    return Th_bps
end

function Performance(L::PHY)
    Th_bps = max_PHY_Th(L)
    #Th = min(L.λ,1/(14*L.T_mu_s))
    Th = L.λ
    Delay = (L.Tproc2/2)+(L.TBS/Th_bps)+(L.radio/299792458)+(L.Tproc1/2)
    return Th, Delay
end

mutable struct Source
    λ_raw # Traffic that arrives to the source (It could be bursty)
    λ_source # Traffic given from the source
    B # Buffer Size
    traffic_shaping_enable # Flag that indicates that the source smooths the traffic
    function Source(input_rate, traffic_shapped = 256, B = 50, flag=false)
        return new(input_rate, traffic_shapped, B, flag)
    end 
end 

function Performance(S::Source)
    if S.traffic_shaping_enable
        # Model
        ρ = S.λ_raw/S.λ_source
        N = ρ*(1-((S.B+1)*ρ^(S.B))+(S.B*ρ^(S.B+1)))/((1-ρ^(S.B+1))*(1-ρ)) # Queue Size
        Th = min(S.λ_raw,S.λ_source)
        W = N/Th

        return Th, W
    else
        return S.λ_raw, 0
    end
end

mutable struct CrossLayer
    S::Source
    L2_4::SDAP
    L2_3::PDCP
    L2_2::RLC_AM
    L2_1::MAC
    L1::PHY
    CL_TH
    CL_Delay
    function CrossLayer(S::Source,L2_4::SDAP,L2_3::PDCP,L2_2::RLC_AM,L2_1::MAC,L1::PHY)
        max_SDAP = L2_4.μ_SDAP
        max_PDPC = (L2_3.μ1_PDCP*L2_3.μ2_PDCP*L2_3.μ3_PDCP)/((L2_3.μ2_PDCP*L2_3.μ3_PDCP)+(L2_3.μ1_PDCP*L2_3.μ3_PDCP)+(L2_3.μ1_PDCP*L2_3.μ2_PDCP))
        max_RLC = min(L2_2.λ_header,L2_2.λ_ACK)
        #max_MAC = min(L2_1.μ_MAC,(1/L2_1.TTI))
        max_MAC = min(L2_1.μ_MAC,L2_1.μ_TX_MAC)
        S.λ_source = floor(min(max_SDAP,max_PDPC,max_RLC,max_MAC))
        tmp = new(S,L2_4,L2_3,L2_2,L2_1,L1,0,0)
        tmp_performance = Performance(tmp)
        tmp.CL_TH = tmp_performance[1]
        tmp.CL_Delay = tmp_performance[2]
        return tmp
    end
end

function Performance(cl::CrossLayer)
    # Function for evaluating the Performance of the CrossLayer design
    s_performance = Performance(cl.S)
    sdpa_performance = Performance(cl.L2_4)
    pdpc_performance = Performance(cl.L2_3)
    rlc_performance = Performance(cl.L2_2)
    mac_performance = Performance(cl.L2_1)
    phy_performance = Performance(cl.L1)
    Th = min(s_performance[1],sdpa_performance[1],pdpc_performance[1],rlc_performance[1],mac_performance[1],phy_performance[1])
    Delay = s_performance[2]+sdpa_performance[2]+pdpc_performance[2]+rlc_performance[2]+mac_performance[2]+phy_performance[2]
    return Th, Delay
end
