import pandas as pd


# Step 1: Create the fatality rate dataframe (dt1)
# Define age-specific fatality rates based on age groups
def get_fatality_rate(age):
    if 0 <= age <= 17:
        return 20 / 1_000_000
    elif 18 <= age <= 49:
        return 500 / 1_000_000
    elif 50 <= age <= 64:
        return 6000 / 1_000_000
    elif age >= 65:
        return 90000 / 1_000_000
    else:
        return 0

# Create dt1 with Age 0 to 100
dt1 = pd.DataFrame({'Age': list(range(0, 101))})
dt1['FatalityRate'] = dt1['Age'].apply(get_fatality_rate)

# Step 2: Read in WorldDemographics.csv (dt2)
dt2 = pd.read_csv('/Users/eva/Desktop/Python/WorldDemographics.csv')  # adjust path if needed
dt2.columns = dt2.columns.str.strip()  # Strip whitespace from column names
print("CSV Columns:", dt2.columns.tolist())  # 👈 THIS LINE IS THE IMPORTANT ONE

# Step 3: Join dt1 and dt2 on Age
merged = dt2.merge(dt1, on='Age')



# Step 4: Calculate expected deaths for each age
merged['ExpectedDeaths'] = merged['#Alive'] * merged['FatalityRate']
print(dt2.columns)


# Step 5: Group by country to get total population and expected deaths
dt3 = merged.groupby('PopulationID').agg(
    TotalPopulation=('#Alive', 'sum'),
    TotalExpectedDeaths=('ExpectedDeaths', 'sum')
).reset_index()


# Step 6: Calculate percentage of population that would die
dt3['PercentDied'] = (dt3['TotalExpectedDeaths'] / dt3['TotalPopulation']) * 100


# Output result to a CSV file
dt3.to_csv('/Users/eva/Desktop/Python/CovidFatalityByCountry.csv', index=False)
