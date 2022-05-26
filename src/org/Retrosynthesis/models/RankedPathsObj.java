package org.Retrosynthesis.models;

public class RankedPathsObj {
    private Double theoretical_yield;
    private Double atpm_generation;
    private Double nadh_consumption;
    private Double nadph_consumption;
    private Double fva_dif;
    private String pathway;


    public RankedPathsObj(Double theoretical_yield, Double atpm_generation, Double nadh_consumption, Double nadph_consumption, Double fva_dif, String path){
        this.theoretical_yield = theoretical_yield;
        this.atpm_generation = atpm_generation;
        this.nadh_consumption = nadh_consumption;
        this.nadph_consumption = nadph_consumption;
        this.fva_dif = fva_dif;
        this.pathway = path;
    }

    public Double getTheoretical_yield() {
        return theoretical_yield;
    }

    public Double getAtpm_generation() {
        return atpm_generation;
    }

    public String getPathway() {
        return pathway;
    }

    public Double getNadh_consumption() {
        return nadh_consumption;
    }

    public Double getNadph_consumption() {
        return nadph_consumption;
    }

    public Double getFva_dif() {
        return fva_dif;
    }
}
