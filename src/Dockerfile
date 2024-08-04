FROM continuumio/miniconda3
RUN conda install conda-forge::rdkit


WORKDIR /app

COPY . /app

RUN pip install --no-cache-dir -r requirements.txt


ENV PORT=8000

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
LABEL authors="Nika"
