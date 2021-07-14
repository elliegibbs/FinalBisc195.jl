module Assignment07

export normalizeDNA
export composition 
export gc_content
export complement 
export reverse_complement
export fasta_header
export parse_fasta

# # uncomment the following line if you intend to use BioSequences types
# using BioSequences

"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases).
Returns a String.
"""
function normalizeDNA(seq)
    seq = uppercase(string(seq))
    for base in seq
        # note: `N` indicates an unknown base
        occursin(base, "AGCTN") || error("invalid base $base")
    end
    return seq # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end

# Your code here.
# Don't forget to export your functions!

function composition(sequence)
    sequence = normalizeDNA(sequence) # make uppercase string, check invalid bases
    a = c = g = t = n = 0 # sets all 4 variables to `0`

    for base in sequence
        # add 1 to each base as it occurs
        #println(typeof(base))
        if base == 'A'
            a = a+1
        elseif base == 'G'
            g = g+1
        elseif base == 'T'
            t = t+1
        elseif base == 'C'
            c = c+1
        elseif base == 'N'
            n = n+1
        end
    end
    
    bases = Dict('A' => a,
                'G' => g,
                'T' => t,
                'N' => n,
                'C' => c)
   

    return bases
end

function gc_content(sequence)
    # Be sure to use `basecomposition()` in your answer.
    # Note: Since `basecomposition()` already calls `normalizeDNA`,
    # there's no need to call it here.
        ## throw an error if the string contains anything other than ACGT
    basecounts = composition(sequence)
   

    seqlength = length(sequence)
    gs = basecounts['G']
    cs = basecounts['C']


    return (gs + cs) / seqlength 

end

function complement_char(base)
    complements =  Dict( 'A' => 'T',
                  'T' => 'A',
                  'G' => 'C',
                  'C' => 'G',
                  'N' => 'N')
    return complements[base]
end

function complement(base::String)
    base = normalizeDNA(base)
    outcome = []
    for letter in base
        letter = complement_char(letter)
        push!(outcome, letter)
    end
    outcome = join(outcome)
return outcome
end

function reverse_complement(sequence)
    sequence = normalizeDNA(sequence)
    backwards = reverse(sequence)
    backarray = []
    for i in 1:length(backwards)
       comp = complement_char(backwards[i])
       push!(backarray,comp)
       #@info backarray
    end 
    return join(string.(backarray))
end

function fasta_header(header)
    if startswith(header, '>')
        header = string(header)
        header = chop(header, head = 1, tail = 0)   
    return header 
    else
        error("Invalid header (headers must start with '>')")
    end
end

function parse_fasta(path)
    ## Think through the components you need
    ## Does it make sense to define any containers at the beginning?
    ## How will you loop through the file?
    ## What do you need to get from each line?
    head_array = []
    body_array = []
    another_array = []
    
    for line in eachline(path)
        
        if startswith(line, '>')
            header_vector = fasta_header(line)
            #header_vector = Tuple([header_vector])
            push!(head_array, header_vector)
            another_array_string = join(another_array)
            if !isempty(another_array) 
                push!(body_array, another_array_string)
            end
            another_array = []
        else
            body_line = line
            body_line = normalizeDNA(line)
            push!(another_array, body_line)
        end
        
    end
    another_array_string = join(another_array)
    push!(body_array, another_array_string)
    Tuple([body_array])
    return head_array, body_array
end
end # module Assignment07
