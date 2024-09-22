import boto3

lambda_client = boto3.client('lambda')

response = lambda_client.invoke(
    FunctionName='HelloStudentFunction',
    InvocationType='RequestResponse',
    Payload=b'{"names": ["Nick", "Bob", "Mike"]}'
)

response_payload = response['Payload'].read()
print(response_payload)
