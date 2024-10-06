# Code for defining the packet

mutable struct packet
    type # Used to identify the type of packet
    creation_time # Time when the packet was created
    delay # Delay of the packet
    rlc_rtx # Variable that records the number of retransmission of each packet in the RLC layer
    size # Size of the packet
    tbs # Transport Block Size
    RBs # required Resource Block
end

