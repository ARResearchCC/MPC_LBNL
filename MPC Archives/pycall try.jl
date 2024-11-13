# ENV["PYTHON"] = "C:\\Users\\Fred\\anaconda3\\python.exe"
using Pkg
# Pkg.build("PyCall")
using PyCall
# Add the directory containing pytry.py to the Python path
pushfirst!(PyVector(pyimport("sys")."path"), "C:\\Users\\Fred\\Desktop\\PyFMI")

# Import the script as a module
your_script = pyimport("your_script")

# Call the fmu function from the script
res = your_script.fmu(1)
# Import necessary Python libraries
pd = pyimport("pandas")

# Convert the pandas DataFrame to a Julia DataFrame
function pandas_to_julia(py_df)
    columns = py_df.columns.tolist()
    data_dict = Dict{String, Any}()
    for col in columns
        data_dict[col] = py_df[col].values
    end
    df = DataFrame(data_dict)
    select!(df, :Time, Not(:Time))
    return df
end

# Convert the DataFrame
df = pandas_to_julia(res)



