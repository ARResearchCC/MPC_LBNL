using DataFrames, XLSX, GLM, Statistics

# Read Excel file into a dictionary of DataFrames
file_path = "C://Users//Fred//Desktop//PyFMI//HP_COP_data.xlsx"
# Create a dictionary to hold DataFrames for each sheet
dfs = Dict{String, DataFrame}()
    
# Open the Excel file
XLSX.openxlsx(file_path) do workbook
    # Get the names of all sheets
    sheet_names = XLSX.sheetnames(workbook)
    
    # Loop through each sheet name and read it into a DataFrame
    for sheet_name in sheet_names
        sheet = XLSX.gettable(workbook[sheet_name])  # Read each sheet as a table
        dfs[sheet_name] = DataFrame(sheet)           # Convert to DataFrame
    end
end

# Extract the "Heating Data" DataFrame
COP_data_H = dfs["Heating Data"]  # Replace with your actual sheet name
COP_data_C = dfs["Cooling Data"]  # Replace with your actual sheet name

function fit_cop(COP_data::DataFrame)
    # Ensure required columns exist
    required_cols = ["Ta", "T_HP", "Compressor_Speed", "COP"]
    if !all(required_cols .∈ Ref(names(COP_data)))
        error("Missing required columns in the DataFrame.")
    end

    # Convert columns to Float64 if necessary
    for col in required_cols
        COP_data[!, col] = convert(Vector{Float64}, COP_data[!, col])
    end

    # Create squared and interaction terms
    COP_data.Ta_sq = COP_data.Ta .^ 2
    COP_data.T_HP_sq = COP_data.T_HP .^ 2
    COP_data.Comp_Speed_sq = COP_data.Compressor_Speed .^ 2
    COP_data.Ta_T_HP = COP_data.Ta .* COP_data.T_HP
    COP_data.Ta_Comp_Speed = COP_data.Ta .* COP_data.Compressor_Speed
    COP_data.T_HP_Comp_Speed = COP_data.T_HP .* COP_data.Compressor_Speed

    # Fit polynomial regression model with interaction terms
    model = lm(@formula(COP ~ Ta + T_HP + Compressor_Speed + 
                              Ta_sq + T_HP_sq + Comp_Speed_sq + 
                              Ta_T_HP + Ta_Comp_Speed + T_HP_Comp_Speed), COP_data)

    # Predict and calculate metrics
    COP_data.COP_pred = predict(model)
    r2_score = r2(model)
    mse = mean((COP_data.COP .- COP_data.COP_pred) .^ 2)
    rmse = sqrt(mse)

    # Print results
    println("R² Score: ", round(r2_score, digits=4))
    println("MSE: ", round(mse, digits=4))
    println("RMSE: ", round(rmse, digits=4))
    println("Coefficients: ", coef(model))

    return coef(model)
end

