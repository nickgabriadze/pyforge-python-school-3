from collections import defaultdict
from airflow import DAG
from airflow.providers.amazon.aws.operators.s3 import S3CreateObjectOperator
from airflow.operators.python import PythonOperator
from airflow.providers.postgres.operators.postgres import PostgresOperator
import json
import pendulum
import logging
import pandas as pd

from rdkit import Chem
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.Lipinski import NumHDonors, NumHAcceptors
from rdkit.Chem.MolSurf import TPSA


def analyze_molecule(smiles, analyzed_dict):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    analyzed_dict["smiles_structure"].append(smiles)
    analyzed_dict["molecular_weight"].append(MolWt(mol))
    analyzed_dict["log_p"].append(MolLogP(mol))
    analyzed_dict["tpsa"].append(TPSA(mol))
    analyzed_dict["h_donors"].append(NumHDonors(mol))
    analyzed_dict["h_acceptors"].append(NumHAcceptors(mol))
    analyzed_dict["lipinski_pass"].append(all([
        MolWt(mol) <= 500,
        NumHDonors(mol) <= 5,
        NumHAcceptors(mol) <= 10,
        MolLogP(mol) <= 5
    ]))


def generate_json_for_molecules(**context):
    molecules = context['ti'].xcom_pull(task_ids='execute_sql')
    logging.info(molecules)
    analyzed_dict = defaultdict(list)

    for molecule in molecules:
        # I have two columns in database table(id, smiles)
        analyze_molecule(molecule[1], analyzed_dict)


    analyzed_frame = pd.DataFrame(analyzed_dict)
    xlsx_output = 'analyzed_smiles.xlsx'
    with pd.ExcelWriter(xlsx_output, engine='openpyxl') as writer:
        analyzed_frame.to_excel(writer, index=False)



default_args = {
    'owner': 'airflow',
}

with DAG(
        dag_id='molecule_analysis_to_s3',
        default_args=default_args,
        description='Analyze molecules and upload results to S3',
        schedule_interval='@daily',
        schedule=None,
        start_date=pendulum.today(),
) as dag:
    retrieve_molecules = PostgresOperator(
        task_id='execute_sql',
        sql='SELECT * FROM molecules LIMIT 10',
        postgres_conn_id='molecules_database'

    )

    generate_json_task = PythonOperator(
        task_id='generate_json_for_molecules',
        python_callable=generate_json_for_molecules,
    )

    upload_to_s3 = S3CreateObjectOperator(
        task_id="create_object",
        s3_bucket='hw-bucket-3000',
        s3_key='analyzed_smiles.xlsx',
        data='./analyzed_smiles.xlsx',
        replace=True,
    )

retrieve_molecules >> generate_json_task >> upload_to_s3
