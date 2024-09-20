import requests
from bs4 import BeautifulSoup
import os

# URL of the transient table
base_url = 'https://swift.gsfc.nasa.gov/results/transients/'

# Create a session
session = requests.Session()

# Get the HTML content of the base URL
response = session.get(base_url)
response.raise_for_status()

# Parse the HTML content
soup = BeautifulSoup(response.content, 'html.parser')

# Find all table rows (assuming the identifiers are in a table)
table = soup.find('table')
rows = table.find_all('tr')[1:]  # Skip the header row

# Create a directory to store the FITS files
if not os.path.exists('fits_files'):
    os.makedirs('fits_files')

# Function to download a file from a URL
def download_file(url, directory):
    local_filename = os.path.join(directory, url.split('/')[-1])
    with session.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    print(f"Downloaded {local_filename}")

name_list = []
# Iterate through each row and construct the .lc.fits URL
for row in rows:
    columns = row.find_all('td')
    if columns:
        print(columns[5].text)
        if ('HMXB' in (columns[5].text)) or ('LMXB' in columns[5].text) or ('XRB' in columns[5].text):
            identifier = columns[1].text.strip()  # Assuming the identifier is in the first column
            name_list.append(identifier)
            identifier = identifier.replace(' ', '')
            identifier = identifier.replace('+', 'p')
            print(identifier)
            fits_url_weak = f'https://swift.gsfc.nasa.gov/results/transients/weak/{identifier}.lc.fits'
            fits_url_main = f'https://swift.gsfc.nasa.gov/results/transients/{identifier}.lc.fits'
            try:
                download_file(fits_url_weak, 'fits_files_XRB')
            except requests.HTTPError:
                print('Usings Alt Link')
                download_file(fits_url_main, 'fits_files_XRB')
with open("names_2.txt", "w") as output:
    output.write(str(name_list))
print(name_list)
print("All files processed.")

