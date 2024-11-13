import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

wind_curtailment_data = client.get_dataset(
    dataset="caiso_curtailment_aggregated_hourly",
    start="2016-06-30",
    end="2023-09-28",
    columns=["wind_curtailment_mw"],
    tz="US/Pacific",
)
solar_curtailment_data = client.get_dataset(
    dataset="caiso_curtailment_aggregated_hourly",
    start="2016-06-30",
    end="2023-09-28",
    columns=["solar_curtailment_mw"],
    tz="US/Pacific",
)

total_curtailment_data = client.get_dataset(
    dataset="caiso_curtailment_aggregated_hourly",
    start="2016-06-30",
    end="2023-09-28",
    columns=["total_curtailment_mw"],
    tz="US/Pacific",
)

total_load_data = client.get_dataset(
    dataset="caiso_standardized_hourly",
    start="2019-04-05",
    end="2024-04-06",
    columns=[
        "interval_start_utc",
        "interval_end_utc",
        "load.load",
    ],
    tz="US/Pacific",
)

wind_data = client.get_dataset(
    dataset="caiso_standardized_hourly",
    start="2019-04-05",
    end="2024-04-06",
    columns=[
        "interval_start_utc",
        "interval_end_utc",
        "fuel_mix.wind"],
    tz="US/Pacific",
)

solar_data = client.get_dataset(
    dataset="caiso_standardized_hourly",
    start="2019-04-05",
    end="2024-04-06",
    columns=[
        "interval_start_utc",
        "interval_end_utc",
        "fuel_mix.solar"],
    tz="US/Pacific",
)

# Convert interval_start_local to datetime if not already
curtailment_data['interval_start_local'] = pd.to_datetime(curtailment_data['interval_start_local'])
load_data['interval_start_local'] = pd.to_datetime(load_data['interval_start_local'])
solar_data['interval_start_local'] = pd.to_datetime(solar_data['interval_start_local'])
wind_data['interval_start_local'] = pd.to_datetime(wind_data['interval_start_local'])

# Extract year
curtailment_data['year'] = curtailment_data['interval_start_local'].dt.year
load_data['year'] = load_data['interval_start_local'].dt.year
solar_data['year'] = solar_data['interval_start_local'].dt.year
wind_data['year'] = wind_data['interval_start_local'].dt.year

# Group by year and calculate annual sums
annual_curtailment = curtailment_data.groupby('year').sum().reset_index()
annual_load = load_data.groupby('year').sum().reset_index()
annual_solar_production = solar_data.groupby('year').sum().reset_index()
annual_wind_production = wind_data.groupby('year').sum().reset_index()

# Calculate relative curtailment
annual_curtailment['solar_curtailment_relative'] = (annual_curtailment['solar_curtailment_mw'] / annual_solar_production['solar_production_mw']) * 100
annual_curtailment['wind_curtailment_relative'] = (annual_curtailment['wind_curtailment_mw'] / annual_wind_production['wind_production_mw']) * 100
annual_curtailment['total_curtailment_relative'] = (annual_curtailment['solar_curtailment_mw'] + annual_curtailment['wind_curtailment_mw']) / annual_load['total_load_mw'] * 100

# Plotting
fig, ax = plt.subplots(figsize=(14, 8))

width = 0.3  # width of the bars
years = annual_curtailment['year']

ax.bar(years - width, annual_curtailment['solar_curtailment_relative'], width, label='Solar Curtailment Relative to Solar Production')
ax.bar(years, annual_curtailment['wind_curtailment_relative'], width, label='Wind Curtailment Relative to Wind Production')
ax.bar(years + width, annual_curtailment['total_curtailment_relative'], width, label='Total Curtailment Relative to Total Load')

ax.set_xlabel('Year')
ax.set_ylabel('Curtailment (%)')
ax.set_title('Annual Relative Curtailment')
ax.legend()

plt.show()