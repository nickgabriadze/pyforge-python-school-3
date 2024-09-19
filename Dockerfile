FROM continuumio/miniconda3
RUN conda install conda-forge::rdkit


WORKDIR /app

COPY ./src ./src
COPY ./tests ./tests
COPY requirements.txt ./
COPY .env ./


RUN pip install --no-cache-dir -r requirements.txt

RUN pytest tests/test_main.py

ENV PORT=8000

EXPOSE 8000

CMD ["uvicorn", "src.main:app", "--host", "0.0.0.0", "--port", "8000"]
LABEL authors="Nika"