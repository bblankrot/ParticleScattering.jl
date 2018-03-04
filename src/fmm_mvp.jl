function FMM_mainMVP_pre!(output, beta, scatteringMatrices, φs::Vector{Float64}, ids::Vector{Int64}, P, mFMM, pre_agg, translated_sum)
    #@simd does not have a positive effect in 0.6.0
    #calculate matrix-vector product - devectorized with pre-preconditioning

    @inbounds begin
    A_mul_B!(output, mFMM.Znear, beta)
    G = length(mFMM.groups)
    fill!(pre_agg,0.0)
    #first preagg all
    for ig2 = 1:G
        for is = 1:mFMM.groups[ig2].size
            indices = (mFMM.groups[ig2].point_ids[is]-1)*(2*P+1)
            for ii = 1:2*P+1
                for iq = 1:mFMM.Q
                    pre_agg[iq,ig2] += mFMM.Agg[ig2][iq,(is-1)*(2*P+1) + ii]*beta[indices + ii]
                end
            end
        end
    end

    for ig1 = 1:G
        #translate plane waves from ig2 to ig1
        fill!(translated_sum,0.0)
        for ig2 = 1:G
            if isempty(mFMM.Trans[(ig1-1)*G + ig2])
                continue
            else
                for iQ = 1:mFMM.Q
                    translated_sum[iQ] -= mFMM.Trans[(ig1-1)*G + ig2][iQ]*pre_agg[iQ,ig2]
                end
                #minus sign because the real equation is (I-ST)x=b
            end
        end
        #disaggregate from ig1 center to ig1's scatterers, producing -Tx
        for is = 1:mFMM.groups[ig1].size
            for ip = 1:2*P+1
                disagged = 0.0im
                for iq = 1:mFMM.Q
                    disagged += conj(mFMM.Agg[ig1][iq,(is-1)*(2*P+1) + ip])*translated_sum[iq]
                end
                output[(mFMM.groups[ig1].point_ids[is]-1)*(2*P+1) + ip] += disagged
            end
        end
    end
    #multiply by S to produce -S
    #temp can be moved outward, but for now this preallocation prevents most
    #dynamic mem alloc by this MVP
    temp = Array{Complex{Float64}}(2*P+1)
    for ic = 1:length(ids)
        rng = (ic-1)*(2*P+1) + (1:2*P+1)
        #copy values to avoid dynamic allocation
        # copy!(temp,1:2*P+1,output,rng)
        for ip = 1:2*P+1
            temp[ip] = output[(ic-1)*(2*P+1) + ip]
        end
        v = view(output,rng)
        if φs[ic] == 0.0
            A_mul_B!(v, scatteringMatrices[ids[ic]], temp)
        else
            #rotate without matrix
            rotateMultipole!(temp,-φs[ic],P)
            A_mul_B!(v, scatteringMatrices[ids[ic]], temp)
            rotateMultipole!(v,φs[ic],P)
        end
    end
    #add identity matrix (Ix)
    output .+= beta
    end #inbounds
    return output
end

function FMM_mainMVP_pre2!(output, beta, scatteringMatrices, φs::Vector{Float64}, ids::Vector{Int64}, P, mFMM, pre_agg, translated_sum)
    #calculate matrix-vector product - partially vectorized with pre-preconditioning

    @inbounds begin
    A_mul_B!(output, mFMM.Znear, beta)
    G = length(mFMM.groups)
    fill!(pre_agg,0.0)
    #first preagg all
    for ig2 = 1:G
        for is = 1:mFMM.groups[ig2].size
            indices = (mFMM.groups[ig2].point_ids[is]-1)*(2*P+1)
            for ii = 1:2*P+1
                for iq = 1:mFMM.Q
                    pre_agg[iq,ig2] += mFMM.Agg[ig2][iq,(is-1)*(2*P+1) + ii]*beta[indices + ii]
                end
            end
        end
    end

    for ig1 = 1:G
        #translate plane waves from ig2 to ig1
        fill!(translated_sum,0.0)
        for ig2 = 1:G
            if isempty(mFMM.Trans[(ig1-1)*G+ig2])#FMMnear(ig1,ig2) - either compute distance or check if Trans is empty
                continue
            else
                for iQ = 1:mFMM.Q
                    translated_sum[iQ] -= mFMM.Trans[(ig1-1)*G + ig2][iQ]*pre_agg[iQ,ig2]
                end
                #minus sign because the real equation is (I-ST)x=b
            end
        end
        #disaggregate from ig1 center to ig1's scatterers, producing -Tx
        for is = 1:mFMM.groups[ig1].size
            for ip = 1:2*P+1
                disagged = 0.0im
                for iq = 1:mFMM.Q
                    disagged += conj(mFMM.Agg[ig1][iq,(is-1)*(2*P+1) + ip])*translated_sum[iq]
                end
                output[(mFMM.groups[ig1].point_ids[is]-1)*(2*P+1) + ip] += disagged
            end
        end
    end
    #multiply by S to produce -ST
    #temp can be moved outward, but for now this preallocation prevents most
    #dynamic mem alloc by this MVP
    temp = Array{Complex{Float64}}(2*P+1)
    for ic = 1:length(ids)
        rng = (ic-1)*(2*P+1) + (1:2*P+1)
        #copy values to avoid dynamic allocation
        for ip = 1:2*P+1
            temp[ip] = output[(ic-1)*(2*P+1) + ip]
        end
        v = view(output,rng)
        if φs[ic] == 0.0
            A_mul_B!(v, scatteringMatrices[ids[ic]], temp)
        else
            #rotate without matrix
            rotateMultipole!(temp,-φs[ic],P)
            A_mul_B!(v, scatteringMatrices[ids[ic]], temp)
            rotateMultipole!(v,φs[ic],P)
        end
    end
    #add identity matrix (Ix)
    output .+= beta
    end #inbounds
    return output
end
