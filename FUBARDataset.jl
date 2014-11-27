function loadDataset(gridInfoFile)
    numGridPoints = 0
    numSites = 0

    open(gridInfoFile,"r") do f
        startGrid=false
        startConditionals=false
        startScalers=false
        i=1
        for line in eachline(f)
            if ismatch(r"^{{", line)
                startGrid=true
                i=1
            end
            if ismatch(r"^}", line)
                startGrid=false
            end
            if startGrid
                spl = split(replace(line, r"[{}\s]+", ""), ",")
                numGridPoints += 1
                i += 1
            end

            if ismatch(r"^{\"conditionals\":{", line)
                startConditionals=true
                i=1
                continue
            end
            if ismatch(r"^}", line)
                startConditionals=false
            end
            if startConditionals
                spl = split(replace(line, r"[{}\s]+", ""), ",")
                numSites=length(spl)
            end

            if ismatch(r"^\"scalers\":{", line)
                startScalers=true
                i=1
                continue
            end
        end
    end
    grid = Array(Float64, numGridPoints, 2)
    conditionals = Array(Float64, numGridPoints, numSites)
    scalars = Array(Float64, numSites)

    open(gridInfoFile,"r") do f
        startGrid=false
        startConditionals=false
        startScalers=false
        i=1
        for line in eachline(f)
            if ismatch(r"^{{", line)
                startGrid=true
                i=1
            end
            if ismatch(r"^}", line)
                startGrid=false
            end
            if startGrid
                spl = split(replace(line, r"[{}\s]+", ""), ",")
                for j = 1:length(spl)
                    grid[i,j] = float(spl[j])
                end
                i += 1
            end

            if ismatch(r"^{\"conditionals\":{", line)
                startConditionals=true
                i=1
                continue
            end
            if ismatch(r"^}", line)
                startConditionals=false
            end
            if startConditionals
                spl = split(replace(line, r"[{}\s]+", ""), ",")
                for j = 1:numSites
                    conditionals[i,j] = float(spl[j])
                end
                i += 1
            end

            if ismatch(r"^\"scalers\":{", line)
                startScalers=true
                i=1
                continue
            end
            if ismatch(r"^}", line)
                startScalers=false
            end

            if startScalers
                spl = split(replace(line, r"[{}\s]+", ""), ",")
                for j = 1:numSites
                    scalars[j] = float(spl[j])
                end
                i += 1
            end

        end
    end

    return numGridPoints, numSites, grid, conditionals, scalars
end

