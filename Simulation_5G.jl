# Code used to simulate the 5G NR Protocol stack
include("Packet.jl")
using Distributions
using Match
using MAT
include("Protocols_5G.jl")

# Struct used to definte the radio channel
mutable struct radioChannel
    # The radio channels is described in the Time and Frequency Domains
    channel #
    nResources # Number of available resources
end

mutable struct CrossLayer_Simulation
    next_event # Next Event Identifier
    next_event_time # Time of the next event
    sim_time # Variable used to control simulation time
    sim_limit # Variable that define the simulation time limit
    transmited_packets # Transmited packets array
    arrival_times # Vector that store the packets arrival times
    next_arrival # Time of the next packet arrival

    # The following values are related to the specific 5G NR protocol stack
    stack5GNR::CrossLayer
    # The values associated to the SDAP protocol are
    sdap_buffer # Variable that represent the SDAP buffer
    sdap_busy # Variable that indicated that the SDAP protocol is busy
    next_sdap_pdu # Time when the next SDAP PDU is created

    # Variables related to the PDCP Layer
    pdcp_buffer # Buffer of the PDCP Layer
    pdcp_busy # Flag that indicated that the PDCP layer is processing a packet
    next_pdcp_pdu # Time when the next PDCP PDU is created

    # Variables required for the RLC simulation
    rlc_buffer # Buffer where the RLC SDU are stored
    rlc_hdr_busy  # Flag that indicated that the RLC protocols is creating a header
    next_rlc_pdu # Time when the next RLC PDU will be created
    rlc_header_buffer # Buffer where the RLC PDUs wait for a room in the Transmission Window
    rlc_tx_windows # RLC Transmission Window
    next_rlc_transmission # Time when the next RLC PDU will be transmitted
    rlc_retransmission_buffer # RLC Retransmission Buffer

    # Variables related to the MAC Layer
    mac_buffer # MAC reception buffer
    mac_hdr_busy # Flag that indicates the MAC protocol is creating a header
    next_mac_pdu # Time when the next MAC PDU will be created
    mac_tx_buffer # Buffer where the packets wait for the their transmission
    next_mac_tx # Time when the next MAC tranmission will ocurr
    TTI # Time Transmission Interval
    radio_channel::radioChannel

    # Variables related to the PHY Layer
    phy_buffer # Struct where the packets in the PHY layer are stored

    function CrossLayer_Simulation(next_event,next_event_time,sim_time,sim_limit,transmited_packets,arrival_times,next_arrival,stack5GNR::CrossLayer,sdap_buffer,sdap_busy,next_sdap_pdu,pdcp_buffer,pdcp_busy,next_pdcp_pdu,rlc_buffer,rlc_hdr_busy,next_rlc_pdu,rlc_header_buffer,rlc_tx_windows,next_rlc_transmission,rlc_retransmission_buffer,mac_buffer,mac_hdr_busy,next_mac_pdu,mac_tx_buffer,TTI)
        nsymbol = Int(round((TTI)/(stack5GNR.L1.T_mu_s))) # Number of symbols that compose the TTI
        nRB = stack5GNR.L1.N_BW_j_mu # Number of RB available in the bandwidth [TS 38.101]
        channel = zeros(Int8,nRB,nsymbol) # The channes is defined as a bitmap, that indicates if a resource is free or not
        channel_resources = nsymbol*nRB # Value that indicates the number of resource in the channel during the TTI
        radio_channel = radioChannel(channel,channel_resources)
        tti_times = 0:TTI:sim_limit
        return new(next_event,next_event_time,sim_time,sim_limit,transmited_packets,arrival_times,next_arrival,stack5GNR::CrossLayer,sdap_buffer,sdap_busy,next_sdap_pdu,pdcp_buffer,pdcp_busy,next_pdcp_pdu,rlc_buffer,rlc_hdr_busy,next_rlc_pdu,rlc_header_buffer,rlc_tx_windows,next_rlc_transmission,rlc_retransmission_buffer,mac_buffer,mac_hdr_busy,next_mac_pdu,mac_tx_buffer,Inf,TTI,radio_channel,[])
    end
end

function increase_Delay(cls::CrossLayer_Simulation,δ)
    # Method used for increasing the packet delays
    # Delay for the SDAP Packets
    for sdap_sdu in cls.sdap_buffer
        sdap_sdu.delay += δ
    end
    # Delay for the PDCP Packets
    for pdcp_sdu in cls.pdcp_buffer
        pdcp_sdu.delay += δ
    end
    # Delay for the RLC Packets
    for rlc_sdu in cls.rlc_buffer
        rlc_sdu.delay += δ
    end
    for rlc_pdu in cls.rlc_header_buffer
        rlc_pdu.delay += δ
    end
    for rlc_packet in cls.rlc_tx_windows
        rlc_packet.delay += δ
    end
    for rlc_packet in cls.rlc_retransmission_buffer
        rlc_packet.delay += δ
    end
    # Delay for the MAC Packets
    for mac_packet in cls.mac_buffer
        mac_packet.delay += δ
    end
    for mac_packet in cls.mac_tx_buffer
        mac_packet.delay += δ
    end
end

function arrival(cls::CrossLayer_Simulation)
    # Function that handles the packet arrival to the stack
    # First, the time elapsed since last event is calculated
    δ = cls.next_event_time-cls.sim_time
    # During this time, th packets have been waiting
    increase_Delay(cls,δ)
    # The packet is created 
    pkt = packet("URLLC",cls.next_event_time,Performance(cls.stack5GNR.S)[2],0,100,0,0)
    append!(cls.arrival_times,cls.next_event_time)

    # Simulation clock is updated
    cls.sim_time += δ
    # The time of the next arrival is calculated
    cls.next_arrival = cls.sim_time + (1/Performance(cls.stack5GNR.S)[1]) #rand(Exponential(1/Performance(cls.stack5GNR.S)[1]))

    # Each time ocurrs an arrival, this is added to the SDAP layer
    append!(cls.sdap_buffer,[pkt])

    # Now, we identify if the SDAP is processing a packet or not
    if !(cls.sdap_busy) & (length(cls.sdap_buffer)>0)
        # If the SDAP layer is idle, the processing is started
        cls.next_sdap_pdu = cls.sim_time+rand(Exponential(1/cls.stack5GNR.L2_4.μ_SDAP))
        cls.sdap_busy = true
    end
end

function SDAP_PDU(cls::CrossLayer_Simulation)
    # Method used for generate the SDAP_PDU
    if cls.sdap_busy
        # The delay packet is increased
        δ = cls.next_event_time-cls.sim_time
        increase_Delay(cls,δ)
        # Creation of the SDAP PDU
        sdap_pdu = popfirst!(cls.sdap_buffer)
        sdap_pdu.size += 1 # a byte of header is added
        # Once the PDU is created, it is sended to the PDCP Layer
        append!(cls.pdcp_buffer,[sdap_pdu])
        # The simulation time is updated
        cls.sim_time += δ
        # Now, we verify that the next packet begin its processing
        if length(cls.sdap_buffer) > 0
            cls.next_sdap_pdu = cls.sim_time+rand(Exponential(1/cls.stack5GNR.L2_4.μ_SDAP))
        else
            cls.next_sdap_pdu = Inf
            cls.sdap_busy = false
        end
        # The processing of the PDCP PDU is started
        if (length(cls.pdcp_buffer) > 0) & !(cls.pdcp_busy)
            cls.next_pdcp_pdu = cls.sim_time+(rand(Exponential(1/cls.stack5GNR.L2_3.μ1_PDCP))+rand(Exponential(1/cls.stack5GNR.L2_3.μ2_PDCP))+rand(Exponential(1/cls.stack5GNR.L2_3.μ3_PDCP)))
            cls.pdcp_busy = true
        end
    end
end

function PDCP_PDU(cls::CrossLayer_Simulation)
    # Method used for generate the PDCP_PDU
    if cls.pdcp_busy
        # The packet delay is updated
        δ = cls.next_event_time-cls.sim_time
        increase_Delay(cls,δ)
        # Creation of the PDCP PDU
        pdcp_pdu = popfirst!(cls.pdcp_buffer)
        pdcp_pdu.size += 6 # A 6-byte header is added
        # Once the PDU is created, it is sended to the RLC Layer
        append!(cls.rlc_buffer,[pdcp_pdu])
        # Simulation time is updated
        cls.sim_time += δ
        # Estimation of the next PDU creation
        if length(cls.pdcp_buffer) > 0
            cls.next_pdcp_pdu = cls.sim_time+(rand(Exponential(1/cls.stack5GNR.L2_3.μ1_PDCP))+rand(Exponential(1/cls.stack5GNR.L2_3.μ2_PDCP))+rand(Exponential(1/cls.stack5GNR.L2_3.μ3_PDCP)))
        else
            cls.next_pdcp_pdu = Inf
            cls.pdcp_busy = false
        end
        # The process of the RLC PDU is started
        if (length(cls.rlc_buffer) > 0) & !(cls.rlc_hdr_busy)
            cls.next_rlc_pdu = cls.sim_time+rand(Exponential(1/cls.stack5GNR.L2_2.λ_header))
            cls.rlc_hdr_busy = true
        end
    end
end

function RLC_PDU(cls::CrossLayer_Simulation)
    # Method used for generate the RLC_PDU
    # This process is a store-and-forward
    if cls.rlc_hdr_busy
        # The packet delay is updated
        δ = cls.next_event_time-cls.sim_time
        increase_Delay(cls,δ)
        # Creation of the RLC PDU
        rlc_pdu = popfirst!(cls.rlc_buffer)
        rlc_pdu.size += 2
        # The created RLC PDU is sended to the RLC Header Buffer
        append!(cls.rlc_header_buffer,[rlc_pdu])
        # Simulation clock is updated
        cls.sim_time += δ
        # First, the next RLC PDU creation time is estimated
        if length(cls.rlc_buffer) > 0
            cls.next_rlc_pdu = cls.sim_time + rand(Exponential(1/cls.stack5GNR.L2_2.λ_header))
        else
            cls.next_rlc_pdu = Inf
            cls.rlc_hdr_busy = false
        end
        # Now, the process of RLC is different. The PDU are not directly sended to the following layer
        # They are stored in a transmission buffer, where they wait for transmission opportunity
        # This buffer is the rlc_hdr_buffer, if the number of current transmitted packets is less 
        # than the Transmission Windows Size, it will be added to the TX Windows
        while ((length(cls.rlc_header_buffer)>0) & ((length(cls.rlc_tx_windows)+length(cls.rlc_retransmission_buffer))<cls.stack5GNR.L2_2.windows_size))
            # The packets on the intermediate buffer are sended to the Transmission Window
            push!(cls.rlc_tx_windows,popfirst!(cls.rlc_header_buffer))
        end
        # In case there are no active tranmission, the trasmission is started
        if (cls.next_rlc_transmission==Inf) & (length(cls.rlc_tx_windows)>0)
            cls.next_rlc_transmission = cls.sim_time + rand(Exponential(1/(cls.stack5GNR.L2_2.λ_ACK/cls.stack5GNR.L2_2.p)))
        end
    end
end

function RLC_Transmit(cls::CrossLayer_Simulation)
    # Method used for RLC transmission
    # First, the delay time is updated
    δ = cls.next_event_time-cls.sim_time
    increase_Delay(cls,δ)
    # The packet to be transmitted is removed from the queue
    rlc_pdu = popfirst!(cls.rlc_tx_windows)
    # The packet tranmission counter is updated
    rlc_pdu.rlc_rtx += 1
    # Each packet has a p probability of been transmitted correctly
    if rand() <= cls.stack5GNR.L2_2.p # Correct transmission
        push!(cls.mac_buffer,rlc_pdu)
    else 
        push!(cls.rlc_retransmission_buffer,rlc_pdu)
    end
    # Now, the packets in the retrans buffer and the header buffer are 
    # added to the tranmission window
    if length(cls.rlc_tx_windows) == 0
    # First, the packets of the retransmission buffer are renewed
        while (length(cls.rlc_retransmission_buffer)>0) & ((length(cls.rlc_tx_windows)+length(cls.rlc_retransmission_buffer))<=cls.stack5GNR.L2_2.windows_size)
            rtx_pkt = popfirst!(cls.rlc_retransmission_buffer)
            if rtx_pkt.rlc_rtx <= cls.stack5GNR.L2_2.max_ret
                push!(cls.rlc_tx_windows,rtx_pkt)
            end
        end
    end
    # Follow, the packets in the header buffer are transmitted to the transmission buffer
    while (length(cls.rlc_header_buffer)>0) & ((length(cls.rlc_tx_windows)+length(cls.rlc_retransmission_buffer))<=cls.stack5GNR.L2_2.windows_size)
        push!(cls.rlc_tx_windows,popfirst!(cls.rlc_header_buffer))
    end
    # The next packet is transmitted
    if length(cls.rlc_tx_windows)>0
        cls.next_rlc_transmission = cls.sim_time + rand(Exponential(1/(cls.stack5GNR.L2_2.λ_ACK/cls.stack5GNR.L2_2.p)))
    else
        cls.next_rlc_transmission = Inf
    end
    # The process of the MAC PDU creation is started
    if (length(cls.mac_buffer)>0) & !(cls.mac_hdr_busy)
        cls.next_mac_pdu = cls.sim_time + rand(Exponential(1/cls.stack5GNR.L2_1.μ_MAC))
        cls.mac_hdr_busy = true
    end
end

function MAC_PDU(cls::CrossLayer_Simulation)
    # Method used for generation of the MAC PDU
    if cls.mac_hdr_busy
        # Estimation of time step
        δ = cls.next_event_time-cls.sim_time
        increase_Delay(cls,δ)
        # Creation of the MAC PDU
        mac_pdu = popfirst!(cls.mac_buffer)
        mac_pdu.size += 2
        # The PDU is added to an intermediate buffer
        # The resource required by the packet for transmission are calculated based on the CQI (1:15)
        channel_vars = matread("Packet_Data_arr.mat")
        mac_pdu.tbs = channel_vars["tbs_arr"][cls.stack5GNR.L1.CQI]
        mac_pdu.RBs = channel_vars["RBs_arr"][cls.stack5GNR.L1.CQI]
        append!(cls.mac_tx_buffer,[mac_pdu])
        # The simulation clock is updated
        cls.sim_time += δ
        # Once the packet has been added to the Transmission Buffer, the next MAC transmission time is estimated
        cls.next_mac_tx = ceil(cls.sim_time/cls.TTI)*TTI
        # Now, in case there are more packets waiting in the MAC Buffer
        if length(cls.mac_buffer) >0
            cls.next_mac_pdu = cls.sim_time + rand(Exponential(1/cls.stack5GNR.L2_1.μ_MAC))
        else
            cls.next_mac_pdu = Inf
            cls.mac_hdr_busy = false
        end
    end
end

function MAC_TX(cls::CrossLayer_Simulation)
    defineResources(cls)
    # The delay of the packets are increased
    δ = cls.next_event_time-cls.sim_time
    increase_Delay(cls,δ)
    # Method used for transmitting packets from MAC Layer
    # Those packets are transmitted in the next transmission opportunity
    # This opportunities are defined by the TTI (Transmission Time Interval)
    if length(cls.mac_tx_buffer) > 0
        rem_buff = []
        #There are at least a packet in the Transmission Buffer
        for pkt in cls.mac_tx_buffer
            if cls.radio_channel.nResources >= pkt.RBs
                # There are available resources for transmit the packet
                pkt = assignResources(cls,pkt)
                append!(cls.phy_buffer,[pkt])
                #append!(cls.transmited_packets,[pkt])
                append!(rem_buff,[pkt])
            end
        end
        # Removing the packets that are sended to the PHY Layer
        for pkt in rem_buff
            locations = findall(x->x==pkt,cls.mac_tx_buffer)
            deleteat!(cls.mac_tx_buffer,locations)
        end
        # The transmission in the PHY layer ocurrs rightaway
        PHY_TX(cls)
    end
    # The simulation time is updated
    cls.sim_time += δ
    # The next transmission time is updated
    if length(cls.mac_tx_buffer) > 0
        cls.next_mac_tx = cls.sim_time + cls.TTI # Transmission is attempted in the next TTI
    else 
        cls.next_mac_tx = Inf
    end
end

function defineResources(cls::CrossLayer_Simulation)
    # Code used to define the channel utilization
    channel_size = size(cls.radio_channel.channel)
    # In this case the channel is free
    cls.radio_channel.channel = zeros(Int8,channel_size[1],channel_size[2])
    cls.radio_channel.nResources = channel_size[1]*channel_size[2]
    allocateResources(cls,0.9)
end

function allocateResources(cls::CrossLayer_Simulation,utilization)
    channel_size = size(cls.radio_channel.channel)
    rb_no = channel_size[1]*channel_size[2]
    sample_channel = rand(Binomial(1,utilization),rb_no)
    cls.radio_channel.channel = reshape(sample_channel,channel_size)
    cls.radio_channel.nResources = rb_no-sum(sample_channel)
end

function assignResources(cls::CrossLayer_Simulation,pkt)
    # Method for assigning resources to the packet
    requiredRBs = pkt.RBs
    no_symbol = 0
    # The channel is searched for free resources
    channelSize = size(cls.radio_channel.channel)
    for rb in 1:channelSize[1]
        for symbol in 1:channelSize[2]
            if cls.radio_channel.channel[rb,symbol] == 0
                # This mean that the resource is free
                cls.radio_channel.channel[rb,symbol] = 1 # The resource is marked as occupied
                cls.radio_channel.nResources -= 1
                requiredRBs -= 1
                if requiredRBs == 0
                    no_symbol = symbol
                    break
                end
            end
        end
        if requiredRBs == 0
            break
        end
    end
    # An additional delay is added to the packet
    pkt.delay += (no_symbol-1)*cls.stack5GNR.L1.T_mu_s
    return pkt
end

function PHY_TX(cls::CrossLayer_Simulation)
    # The packets sended to the PHY layer are transmitted rightaway
    phy_c = max_PHY_Th(cls.stack5GNR.L1)
    while length(cls.phy_buffer) > 0
        pkt = popfirst!(cls.phy_buffer)
        pkt.delay += (cls.stack5GNR.L1.Tproc2/2)+(pkt.tbs/phy_c)+(cls.stack5GNR.L1.radio/299792458)+(cls.stack5GNR.L1.Tproc1/2)
        # The sended packets are stored in the transmitted packet struct
        push!(cls.transmited_packets,pkt)
    end
end

function saveResult(cls::CrossLayer_Simulation)
    # Calculate the mean delay of the packets
    delay = 0.0
    for pkt in cls.transmited_packets
        delay += pkt.delay
    end
    delay /= length(cls.transmited_packets)

    th = length(cls.transmited_packets)/cls.sim_time
    matwrite("Simulation_Results"*string(cls.stack5GNR.S.λ_raw)*".mat",Dict("arr_times"=>mean(cls.arrival_times),"pkt_mean_delay"=>delay,"transmitted_packets"=>th))
end

function Simulation_Cycle(cls::CrossLayer_Simulation)
    # The simulation cycle is executed while the time limit is not reached yet
    while cls.sim_time <= cls.sim_limit
        event_array = [cls.next_arrival,cls.next_sdap_pdu,cls.next_pdcp_pdu,cls.next_rlc_pdu,cls.next_rlc_transmission,cls.next_mac_pdu,cls.next_mac_tx]
        cls.next_event = argmin(event_array)
        cls.next_event_time = [cls.next_arrival,cls.next_sdap_pdu,cls.next_pdcp_pdu,cls.next_rlc_pdu,cls.next_rlc_transmission,cls.next_mac_pdu,cls.next_mac_tx][cls.next_event]
        # Based on that election, the event method is choosen
        @match cls.next_event begin
            1=> arrival(cls)
            2=> SDAP_PDU(cls)
            3=> PDCP_PDU(cls)
            4=> RLC_PDU(cls)
            5=> RLC_Transmit(cls)
            6=> MAC_PDU(cls)
            7=> MAC_TX(cls)
            _=> return
        end
    end
    saveResult(cls)
end
