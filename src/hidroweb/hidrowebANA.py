import requests
from requests.exceptions import (
    ConnectionError,
    Timeout,
    HTTPError,
    RequestException
)
import geopandas as gpd
import pandas as pd
import os

def get_access_token(
        ID: str,
        PWD : str,
        url : str = "https://www.ana.gov.br/hidrowebservice/EstacoesTelemetricas/OAUth/v1",
        timeout: int = 20) -> str:
    """
    Get access to the telemetering stations from ANA's SNIRH (or another API, if provided).

    id = identification number, obtained through the registration
    pwd = password, obtained through the registration
    timeout = time in seconds
    """
    headers = {
        "accept": "*/*",
        "Identificador": f"{ID}",
        "Senha": f"{PWD}"
    }

    print('Attempting connection to the API.')
    try:
        response = requests.get(
            url,
            headers=headers,
            timeout=timeout
        )

        # Raises HTTPError for 4xx/5xx responses
        response.raise_for_status()

        # Try parsing JSON
        try:
            data = response.json()
            print(f"API access, code {data['code']}: {data['message']}.")
            return data['items']['tokenautenticacao']
        except ValueError:
            print("Response is not JSON:")
            print(response.text)
            return None

    except Timeout:
        print("The request timed out.")
        return None
    except ConnectionError:
        print("Failed to connect to the server.")
        print(response.json()['message'])
        return None
    except HTTPError as e:
        print(f"HTTP error occurred: {e} (status code: {response.status_code})")
        print(response.json()['message'])
        return None
    except RequestException as e:
        print(f"Request error: {e}")
        print(response.json()['message'])
        return None


def get_subbasins(
        auth: str,
        start: str = None,
        end: str = None,
        timeout: int = 30) -> dict:
    """
    Get the info of the registered basins in the system.

    auth = token obtained through get_telemstat_access()
    start = 'yyyy-mm-dd' - can filter the first date of update for the data
    end = 'yyyy-mm-dd' - can filter the last date of update for the data
    timeout = time in seconds
    """

    url = "https://www.ana.gov.br/hidrowebservice/EstacoesTelemetricas/HidroSubBacia/v1"

    headers = {
        "accept": "*/*",
        "Authorization": f"Bearer {auth}"
    }

    try:
        response = requests.get(
            url,
            headers=headers,
            timeout=timeout
        )

        # Raise exception for HTTP 4xx/5xx
        response.raise_for_status()

        # Attempt to parse JSON
        try:
            data = response.json()
            return data
        except ValueError:
            print("Response is not JSON:")
            print(response.text)
            return None

    except Timeout:
        print("The request timed out.")
        return None
    except ConnectionError:
        print("Failed to connect to the server.")
        print(response.json()['message'])
        return None
    except HTTPError as e:
        print(f"HTTP error occurred: {e} (status code: {response.status_code})")
        print(response.json()['message'])
        return None
    except RequestException as e:
        print(f"Request error: {e}")
        print(response.json()['message'])
        return None


def get_telemstat_subbasin(
        auth: str,
        id_subbasin : str,
        start : str = None,
        end : str = None,
        timeout : str = 60,
        save_to_geojson : bool = True,
        output_dir: str = None,
        file_name : str = 'stations_all') -> gpd.GeoDataFrame:
    """
    Get the data from the telemetering stations in a specific basin.
    Filter the stations of a given subbasin.

    auth = token obtained through get_telemstat_access()
    id_subbasin = can be checked from the results of get_subbasins()
    start = 'yyyy-mm-dd' - can filter the first date of update for the data
    end = 'yyyy-mm-dd' - can filter the last date of update for the data
    timeout = time in seconds
    output_dir: Output directory (defaults to ./outputs)
    file_name: Output filename without extension
    """

    # Setup output directory
    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), 'outputs')
        os.makedirs(output_dir, exist_ok=True)


    # Retrieve sub basins in the system, and their respective basin code
    print(f'Attempting to retrieve data for sub-basin {id_subbasin}.')
    try:
        subbasins = get_subbasins(auth)
        for subb in subbasins['items']:
            if subb['codigosubbacia'] == id_subbasin:
                id_basin = subb['Bacia_Codigo']
                break
    except Exception as e:
        print(f"Error retrieving {id_subbasin} data: {e}")
        return None

    url = "https://www.ana.gov.br/hidrowebservice/EstacoesTelemetricas/HidroInventarioEstacoes/v1"

    headers = {
        "accept": "*/*",
        "Authorization": f"Bearer {auth}"
    }

    params = {
    "Código da Bacia": id_basin,
    "Data Atualização Inicial (yyyy-MM-dd)": start,
    "Data Atualização Final (yyyy-MM-dd)": end
    }

    print('Retrieving inventory of telemetering stations in the target sub-basin.')
    try:
        response = requests.get(
            url,
            headers=headers,
            params=params,
            timeout=timeout
        )

        # Raise exception for HTTP 4xx/5xx
        response.raise_for_status()

        # Attempt to parse JSON
        try:
            data = response.json()

            # Transform the data into a geodataframe
            gdf = gpd.GeoDataFrame(
                data["items"],
                geometry=gpd.points_from_xy(
                    [i["Longitude"] for i in data["items"]],
                    [i["Latitude"] for i in data["items"]],
                ),
                crs="EPSG:4326"
                )

            # Filter to only get river gauging station
            gdf = gdf[gdf['Tipo_Estacao'] == 'Fluviometrica']

            # Filter geodataframe for the subbasin only
            gdf = gdf[gdf['Sub_Bacia_Codigo'] == id_subbasin]

            print(f'{len(gdf)} fluvial stations found.')

            # Save the data into a geojson for easy retrieval later
            if save_to_geojson:
                print('Saving retrieved fluvial stations in geojson file.')
                gdf.to_file(os.path.join(output_dir, f'{file_name}.geojson'), driver="GeoJSON")
            
            return gdf

        except ValueError:
            print("Response is not JSON:")
            print(response.text)
            return None

    except Timeout:
        print("The request timed out.")
        return None
    except ConnectionError:
        print("Failed to connect to the server.")
        print(response.json()['message'])
        return None
    except HTTPError as e:
        print(f"HTTP error occurred: {e} (status code: {response.status_code})")
        print(response.json()['message'])
        return None
    except RequestException as e:
        print(f"Request error: {e}")
        print(response.json()['message'])
        return None
    

def get_telemstat_flowseries(
        auth: str,
        id_station: str,
        start: str,
        end: str,
        driver: str = "DATA_LEITURA",
        timeout: int = 60):
    """
    Get the info of the registered basins in the system .

    auth = token obtained through get_telemstat_access()
    id_station = station id code
    start = 'yyyy-mm-dd' - can filter the first date of update for the data
    end = 'yyyy-mm-dd' - can filter the last date of update for the data
    timeout = time in seconds
    """

    url = "https://www.ana.gov.br/hidrowebservice/EstacoesTelemetricas/HidroSerieVazao/v1"

    headers = {
        "accept": "*/*",
        "Authorization": f"Bearer {auth}"
    }

    params = {
    "Código da Estação": id_station,
    "Tipo Filtro Data": driver,
    "Data Inicial (yyyy-MM-dd)": start,
    "Data Final (yyyy-MM-dd)": end
    }

    try:
        response = requests.get(
            url,
            headers=headers,
            params=params,
            timeout=timeout
        )

        # Raise exception for HTTP 4xx/5xx
        response.raise_for_status()

        # Attempt to parse JSON
        try:
            data = response.json()
            print(f"Status of retrieved data for station {id_station}, code {data['code']}: {data['message']}.")
            return data
        except ValueError:
            print("Response is not JSON:")
            print(response.text)
            return None

    except Timeout:
        print("The request timed out.")
        return None
    except ConnectionError:
        print("Failed to connect to the server.")
        # print(response.json()['message'])
        return None
    except HTTPError as e:
        print(f"HTTP error occurred: {e} (status code: {response.status_code})")
        # print(response.json()['message'])
        return None
    except RequestException as e:
        print(f"Request error: {e}")
        # print(response.json()['message'])
        return None


def get_telemstat_flowseries_all(
        auth : str,
        gdf_subbasin_stations : gpd.GeoDataFrame,
        start : str = '2025-01-01',
        end : str = '2025-12-01'
        ) -> dict:
    """
    Get the flow series data for all the stations in a given basin.
    No more than 365 days between start date and end data should be used.

    auth = token obtained through get_telemstat_access()
    gdf_subbasin_stations = GeoDataFrame with all stations in the sub-basin
    start = 'yyyy-mm-dd' - can filter the first date of last measurement in the station
    end = 'yyyy-mm-dd' - can filter the first date of last measurement in the station
    """

    station_codes = list(gdf_subbasin_stations['codigoestacao'])

    # Loop through all stations and retrieve the data of river flow series for station
    stations_flow_series = {}
    for id in station_codes:
        stations_flow_series[id] = get_telemstat_flowseries(auth, int(id), start, end)

    return stations_flow_series


def filter_and_save_successfull_flowseries(
        dict_retrieved_flowseries : dict,
        gdf_subbasinstations : gpd.GeoDataFrame,
        save_to_geojson : bool = True,
        output_dir: str = None,
        file_name : str = 'stations_flowseries') -> gpd.GeoDataFrame:

    """
    After retrieving all the flow series, filter for only those which actually have data.
    Save the filtered result in a csv.
    Return the geodataframe of this data.

    dict_retrieved_flowseries: output from get_telemstat_flowseries_all()
    gdf_subbasinstations: GeoDataFrame of target sub-basin being analysed
    output_dir: Output directory (defaults to ./outputs)
    file_name: Output filename without extension
    """

    # Setup output directory
    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), 'outputs')
        os.makedirs(output_dir, exist_ok=True)

    gdf_stats = pd.DataFrame()

    for st in dict_retrieved_flowseries:
        # Add a check to ensure the retrieved data for the station is not None

        try:
            if len(dict_retrieved_flowseries[st]['items']) > 0:
                gdf_temp = pd.DataFrame(dict_retrieved_flowseries[st]['items'])
                gdf_stats = pd.concat([gdf_stats, gdf_temp])
        except Exception as e:
            print(f'Error trying to access values of station {st} ({e}).')

    if len(gdf_stats) == 0:
        print("No data retrieved.")
        return None

    # Filter all the stations in the basin to only those which data could be retrieved
    gdf_subbasinstations = gdf_subbasinstations[gdf_subbasinstations["codigoestacao"].isin(gdf_stats["codigoestacao"])].copy()
 
    # Then join the dataframes, so the flow series data is given to the geodataframe with the geometries of the stations
    gdf_joined = gdf_subbasinstations.merge(
        gdf_stats,
        on="codigoestacao",
        how="left"
    )

    print(f'{len(gdf_subbasinstations['codigoestacao'].unique())} stations retrieved with data.')

    # Save the data into a geojson for easy retrieval later
    if save_to_geojson:
        print('Saving fluvial stations flow series in geojson file.')
        gdf_joined.to_file(os.path.join(output_dir, f'{file_name}.geojson'), driver="GeoJSON")

    return gdf_joined