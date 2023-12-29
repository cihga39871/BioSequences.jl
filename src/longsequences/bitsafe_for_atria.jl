
#=
Author:  Jiacheng Chuan
Purpose: Edit for Atria trimming program
=#

@inline function isbitsafe(seq::LongDNA{4})
    unsafe_isbitsafe(seq) && 
    @inbounds(seq.data[end]) == 0x0000000000000000 &&
    @inbounds if length(seq.data) > 1 && seq.len % 16 != 0x0000000000000000
        (seq.data[end-1] >> (seq.len % 16 * 4) == 0x0000000000000000)
    else
        true
    end
end
isbitsafe(any) = any

@inline function unsafe_isbitsafe(seq::LongDNA{4})
    length(seq.data) == cld(seq.len, 16) + 1
end

"""
    unsafe_extra_bits_to_zeros!(seq::LongDNA{4})

Caution: use only in bitsafe seq!
"""
@inline function unsafe_extra_bits_to_zeros!(seq::LongDNA{4})
    if !isempty(seq)
        remain = (seq.len % 16)
        @inbounds if remain != 0
            seq.data[end-1] &= ~(0xffffffffffffffff << (remain * 4))
        end
    end
    @inbounds seq.data[end] = 0x0000000000000000
    return seq
end

"""
    bitsafe!(seq::LongDNA{4})

Resize `seq.data` to allow loading a pointer `Ptr{UInt64}` safely at the end of `seq`.

Caution: bitsafe LongDNA{4} may not be compatible on all BioSequences functions, especially those do in-place replacement.
"""
@inline function bitsafe!(seq::LongDNA{4})
    if !unsafe_isbitsafe(seq)
        resize!(seq.data, cld(seq.len, 16) + 1)
    end
    unsafe_extra_bits_to_zeros!(seq)
end

bitsafe!(any) = any

# longsequences/transformations.jl
# """
#     resize!(seq::LongDNA{4}, size::Int[, force::Bool=false])

# It overrides `resize!` in BioSequences. Resize a biological sequence `seq`, to a given `size`. The underlying data is bitsafe.
# """
# @inline function Base.resize!(seq::LongSequence{A}, size::Int, force::Bool=false) where {A}
#     if size < 0
#         throw(ArgumentError("size must be non-negative"))
#     else
#         if force | (seq_data_len(A, size) > seq_data_len(A, length(seq)))
#             resize!(seq.data, seq_data_len(A, size))
#         end
#         seq.len = size
#         bitsafe!(seq)
#     end
# end

# longsequences/transformations.jl
# function reverse_complement!(seq::LongSequence{<:NucleicAcidAlphabet})
#     pred = x -> complement_bitpar(x, Alphabet(seq))
#     reverse_data!(pred, seq.data, seq_data_len(seq) % UInt, BitsPerSymbol(seq))
#     zero_offset!(seq)
#     bitsafe!(seq)
# end

# longsequences/transformations.jl
# function reverse_complement(seq::LongSequence{<:NucleicAcidAlphabet})
#     cp = typeof(seq)(undef, unsigned(length(seq)))
#     pred = x -> complement_bitpar(x, Alphabet(seq))
#     reverse_data_copy!(pred, cp.data, seq.data, seq_data_len(seq) % UInt, BitsPerSymbol(seq))
#     zero_offset!(cp)
#     bitsafe!(cp)
# end
