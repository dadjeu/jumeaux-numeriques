"""Détection simple des cas critiques : une grandeur comparée à un seuil."""


def verifier(valeur, seuil, message):
    """Renvoie un message d'alerte si la valeur dépasse le seuil, sinon None."""
    if valeur > seuil:
        return f"ALERTE : {message} ({valeur:.3g} > seuil {seuil:.3g})"
    return None
