function MH = hmchop(H)
    MH = H;

    if ~isempty(H.D)
        MH.D = mchop(H.D);
    else
        MH.U1 = mchop(H.U1);
        MH.V1 = mchop(H.V1);
        MH.U2 = mchop(H.U2);
        MH.V2 = mchop(H.V2);
        MH.A11 = hmchop(H.A11);
        MH.A22 = hmchop(H.A22);
    end
end
